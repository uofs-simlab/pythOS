SUNDIALS_INSTALL_DIR=sundials_install/

SUNDIALS_LIB_DIR=lib

all: sundials_wrapper.so

sundials_wrapper.so: sundials_wrapper.c
	gcc -shared -o sundials_wrapper.so sundials_wrapper.c -I$(SUNDIALS_INSTALL_DIR)/include -lsundials_core -L$(SUNDIALS_INSTALL_DIR)/$(SUNDIALS_LIB_DIR)
