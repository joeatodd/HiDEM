#Example of how to run HiDEM w/ valgrind to check for memory issues
mpirun -n 4 valgrind --tool=memcheck --leak-check=yes --track-origins=yes ./HiDEM |& tee LOG_valgrind 
