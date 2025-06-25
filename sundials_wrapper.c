#include <cvode/cvode.h>
#include <arkode/arkode.h>

#ifdef __unix__
#define EXPORT
#elif defined(_WIN32) || defined(WIN32)
#define EXPORT __declspec(dllexport)
#endif

EXPORT SUNContext allocate_context() {
  int err;
  SUNContext ctx;

  err = SUNContext_Create(SUN_COMM_NULL, &ctx);

  return ctx;
}

EXPORT MRIStepInnerStepper* allocate_inner_stepper() {
  MRIStepInnerStepper*  stepper;
  stepper = (MRIStepInnerStepper *) malloc(sizeof(MRIStepInnerStepper));

  return stepper;
}
