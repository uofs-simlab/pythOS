SUNDIALS_INCLUDE_DIR=../../sundials/instdir/include/

all: sundials_wrapper.so

sundials_wrapper.so: sundials_wrapper.c
	gcc -shared -o sundials_wrapper.so sundials_wrapper.c -I$(SUNDIALS_INCLUDE_DIR)
