#include "stdlib.h"
#include "stdio.h"
#include "evio.h"
#include "string.h"

#define BUFFER_LENGTH 100000

int main(int argc, char *argv[])
{
    char *ptr;
    char outputf[256], inputf[256];
    int i;
    for(i = 1; i < argc; ++i)
    {
        ptr = argv[i];
        if(*(ptr++) == '-') {
            switch(*(ptr++))
            {
            case 'o':
                strcpy(outputf, argv[++i]);
                break;
            case 'i':
                strcpy(inputf, argv[++i]);
                break;
            default:
                printf("Unkown option!\n");
                exit(1);
            }
        }
    }

    int outHandle, status;
    uint32_t buffer[BUFFER_LENGTH];
    uint32_t event_length;

    FILE *infile = fopen(inputf, "r");
    if(infile == NULL) {
        printf("Cannot find input file \"%s\"", inputf);
        exit(1);
    }

    status = evOpen(outputf, "w", &outHandle);
    if(status != S_SUCCESS) {
        printf("Error in open output file!\n");
        exit(1);
    }

    printf("Start to convert file %s to %s...\n", inputf, outputf);

    int index = 0;
    while(index < BUFFER_LENGTH) {
        status = fscanf(infile, "0x%x\n", &buffer[index]);
        printf("0x%x\n", buffer[index]);
        if(status != 1) break;
        ++index;
    }

    printf("event length: %d\n", buffer[0]);

    evWrite(outHandle, buffer);

    fclose(infile);
    evClose(outHandle);

    return 0;
}
