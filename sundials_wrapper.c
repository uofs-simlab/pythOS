#include <cvode/cvode.h>
#include <arkode/arkode.h>

__declspec(dllexport) SUNContext* allocate_context() {
  SUNContext* ctx;
  ctx = (SUNContext *) malloc(sizeof(SUNContext));
  
  return ctx;
}

__declspec(dllexport) MRIStepInnerStepper* allocate_inner_stepper() {
  MRIStepInnerStepper*  stepper;
  stepper = (MRIStepInnerStepper *) malloc(sizeof(MRIStepInnerStepper));

  return stepper;
}
