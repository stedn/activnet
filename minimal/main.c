/*
** main.c
*/
#include <stdio.h>
#include <stdlib.h>
#include "activnet_gen.h"

int main( int argc, char *argv[] )
{
    if(argc==18){
        
        activnet_gen(atof(argv[1]),atof(argv[2]),atof(argv[3]),atof(argv[4]),atof(argv[5]),atof(argv[6]),
                atof(argv[7]),atof(argv[8]),atof(argv[9]),atof(argv[10]),atof(argv[11]),atof(argv[12]),
                atof(argv[13]),atof(argv[14]),atof(argv[15]),atof(argv[16]),atof(argv[17]));
        
        
    }else{
        printf("Usage: activnet_gen zet L mu kap lc del ups phi psi r sig D Df ls lf tinc tfin \n");
        exit(1);
    }
    return 0;
}
