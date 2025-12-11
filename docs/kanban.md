## TODOs

### Core Development

- ros2 app: do fast imu prop so state is not delayed.

- Visual Odometry (simplified groves/px4 style)

- Multiple hypotheses 

- Particle filter

- real time whiteness check

- better eval: display NEES, NIS, MAE, RMSE

- monte carlo bash script + eval script

- in-house trajectory generator

- no dynamic memory allocation

- Joseph update

- templating: add double and float aliases, rewrite apps to use those. no weird compile time type

### Support & Infrastructure

- single app for closed loop / open loop, with yaml opencv parser

- profile code with valgrind

- write unit tests + test bash script

- Remove dependencies: place eigen in include dir as submodule

- clang tidy: cleanup code

- doxygen: add more latex to documentation, complete docs stubs

- why tests fail with float build? But apps seem to work.

- document apps, add jsbsim version tag

- improve readme, add images & gifs, better plots

# DONE

- generic core templated kalman update func

- template double float filter funcs and classes

- modular update functions, no duplications

- Git hooks

- no matrix inversions!