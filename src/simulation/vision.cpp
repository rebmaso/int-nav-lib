#include "simulation.h"

namespace intnavlib {

#ifdef BUILD_VISION

template<typename T>
Simulation<T>::ViewerInterface::ViewerInterface(const std::string& earth_file_path,
                                                double fx, 
                                                double fy, 
                                                double cx, 
                                                double cy, 
                                                unsigned int width, 
                                                unsigned int height, 
                                                double near, 
                                                double far, 
                                                unsigned int delay_frames,
                                                unsigned int init_delay_frames) :

                                                fx(fx),
                                                fy(fy),
                                                cx(cx),
                                                cy(cy),
                                                width(width),
                                                height(height),
                                                near(near),
                                                far(far),
                                                delay_frames(delay_frames),
                                                init_delay_frames(init_delay_frames) {

        osgEarth::initialize();

        std::vector<std::string> fake_argv_str;
        fake_argv_str.push_back("nav_sim_ol"); // dummy program name
        fake_argv_str.push_back("--nologdepth");
        fake_argv_str.push_back(earth_file_path);

        std::vector<char*> fake_argv;
        for (auto& s : fake_argv_str) {
            fake_argv.push_back(&s[0]);
        }
        int fake_argc = fake_argv.size();
        arguments = std::make_unique<osg::ArgumentParser>(&fake_argc, fake_argv.data());
        
        viewer = new osgViewer::Viewer(*arguments);
        
        // // EGL graphics context for headless rendering
        // // When headless, only works on NVIDIA GPUs it seems
        // // If i call eglinitialize (eg with eglinfo) it fails when headless
        // // TODO: try with the exact same impl from the nvidia example
        // // https://community.khronos.org/t/offscreen-rendering-using-opengl-and-egl/109487/19
        // // https://developer.nvidia.com/blog/egl-eye-opengl-visualization-without-x-server/
        // // https://github.com/omaralvarez/osgEGL/tree/master
        gc = new EGLGraphicsWindowEmbedded(static_cast<int>(width), static_cast<int>(height));
        viewer->getCamera()->setGraphicsContext(gc);
        viewer->getCamera()->setViewport(new osg::Viewport(0, 0, width, height));
        viewer->realize();

        // Setup hidden viewer using pbuffer
        // Should work without a display on intel GPUs
        // https://osg-users.openscenegraph.narkive.com/pIDGiyl0/hidden-viewer
        // traits = new osg::GraphicsContext::Traits;
        // traits->x = 0;
        // traits->y = 0;
        // traits->width = width;
        // traits->height = height;
        // traits->windowDecoration = false;
        // traits->doubleBuffer = false;
        // traits->sharedContext = 0;
        // traits->pbuffer = true; // Set pixel buffer
        // traits->readDISPLAY();
        // gc = osg::GraphicsContext::createGraphicsContext(traits.get());
        // gc->realize();
        // gc->makeCurrent();
        // viewer->getCamera()->setGraphicsContext(gc);
        // viewer->getCamera()->setViewport(new osg::Viewport(0, 0, width, height));

        manip = new SimpleManipulator();
        viewer->setCameraManipulator(manip.get());

        viewer->getCamera()->setSmallFeatureCullingPixelSize(-1.0f);
        viewer->getCamera()->setProjectionResizePolicy(osg::Camera::FIXED); 
        viewer->getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);

        auto map_node = MapNodeHelper().load(*arguments, viewer.get());
        if (map_node.valid() && MapNode::get(map_node))
            viewer->setSceneData(map_node);

        setProj(fx, fy, cx, cy, width, height, near, far);

        // LOG(INFO) << "" << std::endl;
        // LOG(INFO) << "" << "OpenGL Vendor: " << glGetString(GL_VENDOR) << std::endl;
        // LOG(INFO) << "" << "OpenGL Renderer: " << glGetString(GL_RENDERER) << std::endl;
        // LOG(INFO) << "" << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
        // LOG(INFO) << "" << "OpenGL Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
        // LOG(INFO) << "" << std::endl;

        // Allocate and attach a persistent color image.
        color_image = new osg::Image();
        color_image->allocateImage(static_cast<int>(width), static_cast<int>(height), 1, GL_RGB, GL_UNSIGNED_BYTE);
        viewer->getCamera()->attach(osg::Camera::COLOR_BUFFER, color_image.get());

        // Allocate and attach a persistent depth image.
        depth_image = new osg::Image();
        depth_image->allocateImage(static_cast<int>(width), static_cast<int>(height), 1, GL_DEPTH_COMPONENT, GL_FLOAT);
        viewer->getCamera()->attach(osg::Camera::DEPTH_BUFFER, depth_image.get());

}

template<typename T>
void Simulation<T>::ViewerInterface::setPose(const Eigen::Matrix3d &rotation, const Eigen::Vector3d &translation) {

	// Construct 4x4 Eigen transformation matrix
	Eigen::Matrix4d eigenViewMatrix = Eigen::Matrix4d::Identity();

	// Invert here
	eigenViewMatrix.block<3, 3>(0, 0) = rotation.transpose();
	eigenViewMatrix.block<3, 1>(0, 3) = - rotation.transpose() * translation;


	// Flip Y and Z axes
	eigenViewMatrix.row(1) *= -1;
	eigenViewMatrix.row(2) *= -1;

	// Transpose the matrix
	eigenViewMatrix.transposeInPlace();

	// Construct osg::Matrixd 
	osg::Matrixd viewMatrix(
	    eigenViewMatrix(0, 0), eigenViewMatrix(0, 1), eigenViewMatrix(0, 2), eigenViewMatrix(0, 3),
	    eigenViewMatrix(1, 0), eigenViewMatrix(1, 1), eigenViewMatrix(1, 2), eigenViewMatrix(1, 3),
	    eigenViewMatrix(2, 0), eigenViewMatrix(2, 1), eigenViewMatrix(2, 2), eigenViewMatrix(2, 3),
	    eigenViewMatrix(3, 0), eigenViewMatrix(3, 1), eigenViewMatrix(3, 2), eigenViewMatrix(3, 3));
	    
	// Apply the matrix
    	manip->setByInverseMatrix(viewMatrix);
};

template<typename T>
void Simulation<T>::ViewerInterface::setProj(const double & fx, 
	const double & fy, 
	const double & cx, 
	const double & cy, 
	const double & w, 
	const double & h, 
	const double & near, 
	const double & far) {
	osg::Matrixd projMatrix(2 * fx / w, 0.0, (w - 2 * cx) / w, 0.0,
	    0.0, 2 * fy / h, (h - 2 * cy) / h, 0.0,
	    0.0, 0.0, (-far - near) / (far - near), -2.0 * far * near / (far - near),
	    0.0, 0.0, -1.0, 0.0);
	projMatrix.transpose(projMatrix);

	viewer->getCamera()->setProjectionMatrix(projMatrix);
}

template<typename T>
void Simulation<T>::ViewerInterface::startWithInit(const Eigen::Matrix3d &rotation, 
                                        const Eigen::Vector3d &translation) {

	setPose(rotation, translation);
	// LOG(INFO) << "Scene Loading...";
	for(int i = 0; i < init_delay_frames; i++) viewer->frame();
	// LOG(INFO) << "Ready!";                      
	renderThread = std::thread(&ViewerInterface::renderLoop, this);
}

template<typename T>
cv::Mat Simulation<T>::ViewerInterface::getSnapShot() {

	// Wait some frames - producer-consumer model
	std::unique_lock<std::mutex> lk(viewer_mutex);
	viewer_flag = true;
	viewer_cv.wait(lk,[&]{return !viewer_flag;});

	// Read the attached color image, flip it vertically, and encode.
	cv::Mat raw_image;
	raw_image = cv::Mat(color_image->t(), color_image->s(), CV_8UC3, color_image->data()).clone();

	cv::Mat final_image;
	cv::flip(raw_image, final_image, 0);
	cv::cvtColor(final_image, final_image, cv::COLOR_BGR2GRAY);

	return final_image.clone();
}

template<typename T>
void Simulation<T>::ViewerInterface::renderLoop () {
	// Wait for first pose to be set, before looping
	while(!viewer->done()) {
	    viewer->frame();
	    if(viewer_flag) {
		for(size_t i = 0; i<delay_frames; i++) viewer->frame();
		viewer_flag = false;
		viewer_cv.notify_one();
	    }
	}
}

#endif

template struct Simulation<double>;
template struct Simulation<float>;

};