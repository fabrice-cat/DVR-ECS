#include<stdio.h>
#include<stdlib.h>  


void print_symbol_(void); // Be very carreful here, fucntion has be named with _ in the c code while it used without it in fortran
void print_val_(int *);



void print_symbol_(void)
{
   printf("/*"); fflush(stdout);
}

void print_val_(int *i) // Be very carreful here, arguments has to be passed not in value!
{
   printf("\b\b\b\b\b%5d",*i); fflush(stdout);
}
