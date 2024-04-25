#include <cvode/cvode.h>
#include <arkode/arkode.h>

SUNContext* allocate_context() {
  SUNContext* ctx;
  ctx = (SUNContext *) malloc(sizeof(SUNContext));
  
  return ctx;
}

/*void free_cvode(void * cvode_mem) {
  CVodeFree(&cvode_mem);
  }*/

MRIStepInnerStepper* allocate_inner_stepper() {
  MRIStepInnerStepper*  stepper;
  stepper = (MRIStepInnerStepper *) malloc(sizeof(MRIStepInnerStepper));

  return stepper;
}
