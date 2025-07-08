SUNDIALS_INCLUDE_DIR=~/simlab/sundials_install/include/

all: sundials_wrapper.so

sundials_wrapper.so: sundials_wrapper.c
	gcc -shared -o sundials_wrapper.so sundials_wrapper.c -I$(SUNDIALS_INCLUDE_DIR) -lsundials_core -L$(SUNDIALS_INCLUDE_DIR)/../lib
