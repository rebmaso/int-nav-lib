#ifndef EGL_GRAPHICS_WINDOW_EMBEDDED_H
#define EGL_GRAPHICS_WINDOW_EMBEDDED_H

#include <EGL/egl.h>
#include <stdexcept> // For runtime_error
#include <iostream>  // For error messages

// Taken from https://github.com/omaralvarez/osgEGL/tree/master?tab=GPL-3.0-1-ov-file#readme
// See annexed license 
class EGLGraphicsWindowEmbedded : public osgViewer::GraphicsWindowEmbedded
{
protected:
    EGLDisplay  eglDpy;
    EGLint      major, minor;
    EGLint      numConfigs;
    EGLConfig   eglCfg;
    EGLSurface  eglSurf;
    EGLContext  eglCtx;

public:

    EGLGraphicsWindowEmbedded(const int & w, const int & h)
    {
        const EGLint configAttribs[] = {
            EGL_SURFACE_TYPE, EGL_PBUFFER_BIT,
            EGL_BLUE_SIZE, 8,
            EGL_GREEN_SIZE, 8,
            EGL_RED_SIZE, 8,
            EGL_DEPTH_SIZE, 8,
            EGL_RENDERABLE_TYPE, EGL_OPENGL_BIT,
            EGL_NONE
        };

        static const int pbufferWidth = w;
        static const int pbufferHeight = h;

        const EGLint pbufferAttribs[] = {
            EGL_WIDTH, pbufferWidth,
            EGL_HEIGHT, pbufferHeight,
            EGL_NONE,
        };

        // 1. Initialize EGL
        eglDpy = eglGetDisplay(EGL_DEFAULT_DISPLAY);
        if (eglDpy == EGL_NO_DISPLAY) throw std::runtime_error("Failed to get EGL display");

        if (!eglInitialize(eglDpy, &major, &minor)) throw std::runtime_error("Failed to initialize EGL");

        // 2. Select an appropriate configuration
        if (!eglChooseConfig(eglDpy, configAttribs, &eglCfg, 1, &numConfigs))
            throw std::runtime_error("Failed to choose EGL config");

        // 3. Create a surface
        eglSurf = eglCreatePbufferSurface(eglDpy, eglCfg, pbufferAttribs);
        if (eglSurf == EGL_NO_SURFACE) throw std::runtime_error("Failed to create EGL surface");

        // 4. Bind the API
        if (!eglBindAPI(EGL_OPENGL_API)) throw std::runtime_error("Failed to bind EGL API");

        // 5. Create a context and make it current
        eglCtx = eglCreateContext(eglDpy, eglCfg, EGL_NO_CONTEXT, NULL);
        if (eglCtx == EGL_NO_CONTEXT) throw std::runtime_error("Failed to create EGL context");

        // Somehow breaks it
        // if (!eglMakeCurrent(eglDpy, eglSurf, eglSurf, eglCtx)) throw std::runtime_error("Failed to make EGL context current");

        init();
    }

    virtual ~EGLGraphicsWindowEmbedded()
    {
        // 6. Terminate EGL when finished
        if (!eglTerminate(eglDpy)) {
            std::cerr << "Failed to terminate EGL display" << std::endl;
        }
    }

    virtual bool isSameKindAs(const Object* object) const { return dynamic_cast<const EGLGraphicsWindowEmbedded*>(object) != 0; }
    virtual const char* libraryName() const { return "osgViewer"; }
    virtual const char* className() const { return "EGLGraphicsWindowEmbedded"; }

    void init()
    {
        if (valid())
        {
            setState(new osg::State);
            getState()->setGraphicsContext(this);
            getState()->setContextID(osg::GraphicsContext::createNewContextID());
        }
    }

    virtual bool valid() const { return true; }
    virtual bool realizeImplementation() { return true; }
    virtual bool isRealizedImplementation() const { return true; }
    virtual void closeImplementation() {}

    virtual bool makeCurrentImplementation() {
        if (!eglMakeCurrent(eglDpy, eglSurf, eglSurf, eglCtx)) {
            // throw std::runtime_error("Failed to make EGL context current in makeCurrentImplementation");
        }
        return true;
    }

    virtual bool releaseContextImplementation() { return true; }
    virtual void swapBuffersImplementation() {}
    virtual void grabFocus() {}
    virtual void grabFocusIfPointerInWindow() {}
    virtual void raiseWindow() {}
};

#endif
