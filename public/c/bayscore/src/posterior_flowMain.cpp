#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include "fasta-io.h"
#include "zutil.h"

#include "posterior_flow.h"

int main(int argc, char **argv)
{
    return rescorer_main(argc, argv);
}
