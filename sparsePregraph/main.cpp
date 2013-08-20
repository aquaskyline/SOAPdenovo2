#include <stdio.h>
#include <stdlib.h>

extern "C" int call_pregraph_sparse ( int argc, char ** argv );

int main ( int argc, char ** argv )
{
        fprintf ( stderr, "\nVersion 1.0.3: released on July 13th, 2012\nCompile %s\t%s\n\n", __DATE__, __TIME__ );
        call_pregraph_sparse ( argc, argv );
}

