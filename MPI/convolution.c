//
//  convolution.c
//
//
//  Created by Josep Lluis Lerida on 11/03/15.
//
// This program calculates the convolution for PPM images.
// The program accepts an PPM image file, a text definition of the kernel matrix and the PPM file for storing the convolution results.
// The program allows to define image partitions for processing large images (>500MB)
// The 2D image is represented by 1D vector for chanel R, G and B. The convolution is applied to each chanel separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>

// Structure to store image.
struct imagenppm{
    int altura;
    int ancho;
    char *comentario;
    int maxcolor;
    int P;
    int *R;
    int *G;
    int *B;
};
typedef struct imagenppm* ImagenData;

// Structure to store the kernel.
struct structkernel{
    int kernelX;
    int kernelY;
    float *vkern;
};
typedef struct structkernel* kernelData;

//Functions Definition
ImagenData initimage(char* nombre, FILE **fp, int partitions, int halo);
ImagenData duplicateImageData(ImagenData src, int partitions, int halo);

int readImage(ImagenData Img, FILE **fp, int dim, int halosize, long int *position);
int duplicateImageChunk(ImagenData src, ImagenData dst, int dim);
int initfilestore(ImagenData img, FILE **fp, char* nombre, long *position);
int savingChunk(ImagenData img, FILE **fp, int dim, int offset);
int convolve2D(int* inbuf, int* outbuf, int sizeX, int sizeY, float* kernel, int ksizeX, int ksizeY);
void freeImagestructure(ImagenData *src);

//Open Image file and image struct initialization
ImagenData initimage(char* nombre, FILE **fp,int partitions, int halo){
    char c;
    char comentario[300];
    int i=0,chunk=0;
    ImagenData img=NULL;

    /*Opening ppm*/

    if ((*fp=fopen(nombre,"r"))==NULL){
        perror("Error: ");
    }
    else{
        //Memory allocation
        img=(ImagenData) malloc(sizeof(struct imagenppm));

        //Reading the first line: Magical Number "P3"
        fscanf(*fp,"%c%d ",&c,&(img->P));

        //Reading the image comment
        while((c=fgetc(*fp))!= '\n'){comentario[i]=c;i++;}
        for(;i<299;i++)comentario[i]=' ';
        comentario[299]='\0';
        //Allocating information for the image comment
        img->comentario = calloc(strlen(comentario),sizeof(char));
        strcpy(img->comentario,comentario);
        //Reading image dimensions and color resolution
        fscanf(*fp,"%d %d %d",&img->ancho,&img->altura,&img->maxcolor);
        chunk = img->ancho*img->altura / partitions;
        //We need to read an extra row.
        chunk = chunk + img->ancho * halo;
        if ((img->R=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
        if ((img->G=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
        if ((img->B=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
    }
    return img;
}

//Duplicate the Image struct for the resulting image
ImagenData duplicateImageData(ImagenData src, int partitions, int halo){
    int chunk=0;
    //Struct memory allocation
    ImagenData dst=(ImagenData) malloc(sizeof(struct imagenppm));

    //Copying the magic number
    dst->P=src->P;
    //Copying the string comment
    dst->comentario = calloc(strlen(src->comentario),sizeof(char));
    strcpy(dst->comentario,src->comentario);
    //Copying image dimensions and color resolution
    dst->ancho=src->ancho;
    dst->altura=src->altura;
    dst->maxcolor=src->maxcolor;
    chunk = dst->ancho*dst->altura / partitions;
    //We need to read an extra row.
    chunk = chunk + src->ancho * halo;
    if ((dst->R=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
    if ((dst->G=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
    if ((dst->B=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
    return dst;
}

//Read the corresponding chunk from the source Image
int readImage(ImagenData img, FILE **fp, int dim, int halosize, long *position){
    int i=0, k=0,haloposition=0;
    if (fseek(*fp,*position,SEEK_SET))
        perror("Error: ");
    haloposition = dim-(img->ancho*halosize*2);
    for(i=0;i<dim;i++) {
        // When start reading the halo store the position in the image file
        if (halosize != 0 && i == haloposition) *position=ftell(*fp);
        fscanf(*fp,"%d %d %d ",&img->R[i],&img->G[i],&img->B[i]);
        k++;
    }
//    printf ("Readed = %d pixels, posicio=%lu\n",k,*position);
    return 0;
}

//Duplication of the  just readed source chunk to the destiny image struct chunk
int duplicateImageChunk(ImagenData src, ImagenData dst, int dim){
    int i=0;

    for(i=0;i<dim;i++){
        dst->R[i] = src->R[i];
        dst->G[i] = src->G[i];
        dst->B[i] = src->B[i];
    }
//    printf ("Duplicated = %d pixels\n",i);
    return 0;
}

// Open kernel file and reading kernel matrix. The kernel matrix 2D is stored in 1D format.
kernelData leerKernel(char* nombre){
    FILE *fp;
    int i=0;
    kernelData kern=NULL;

    /*Opening the kernel file*/
    fp=fopen(nombre,"r");
    if(!fp){
        perror("Error: ");
    }
    else{
        //Memory allocation
        kern=(kernelData) malloc(sizeof(struct structkernel));

        //Reading kernel matrix dimensions
        fscanf(fp,"%d,%d,", &kern->kernelX, &kern->kernelY);
        kern->vkern = (float *)malloc(kern->kernelX*kern->kernelY*sizeof(float));

        // Reading kernel matrix values
        for (i=0;i<(kern->kernelX*kern->kernelY)-1;i++){
            fscanf(fp,"%f,",&kern->vkern[i]);
        }
        fscanf(fp,"%f",&kern->vkern[i]);
        fclose(fp);
    }
    return kern;
}

// Open the image file with the convolution results
int initfilestore(ImagenData img, FILE **fp, char* nombre, long *position){
    /*Se crea el fichero con la imagen resultante*/
    if ( (*fp=fopen(nombre,"w")) == NULL ){
        perror("Error: ");
        return -1;
    }
    /*Writing Image Header*/
    fprintf(*fp,"P%d\n%s\n%d %d\n%d\n",img->P,img->comentario,img->ancho,img->altura,img->maxcolor);
    *position = ftell(*fp);
    return 0;
}

// Writing the image partition to the resulting file. dim is the exact size to write. offset is the displacement for avoid halos.
int savingChunk(ImagenData img, FILE **fp, int dim, int offset){
    int i,k=0;
    //Writing image partition
    for(i=offset;i<dim+offset;i++){
        fprintf(*fp,"%d %d %d ",img->R[i],img->G[i],img->B[i]);
        k++;
    }
    return 0;
}

// This function free the space allocated for the image structure.
void freeImagestructure(ImagenData *src){

    free((*src)->comentario);
    free((*src)->R);
    free((*src)->G);
    free((*src)->B);

    free(*src);
}

///////////////////////////////////////////////////////////////////////////////
// 2D convolution
// 2D data are usually stored in computer memory as contiguous 1D array.
// So, we are using 1D array for 2D data.
// 2D convolution assumes the kernel is center originated, which means, if
// kernel size 3 then, k[-1], k[0], k[1]. The middle of index is always 0.
// The following programming logics are somewhat complicated because of using
// pointer indexing in order to minimize the number of multiplications.
//
//
// signed integer (32bit) version:
///////////////////////////////////////////////////////////////////////////////
int convolve2D(int* in, int* out, int dataSizeX, int dataSizeY,
               float* kernel, int kernelSizeX, int kernelSizeY)
{
    int i, j, m, n;
    int *inPtr, *inPtr2, *outPtr;
    float *kPtr;
    int kCenterX, kCenterY;
    int rowMin, rowMax;                             // to check boundary of input array
    int colMin, colMax;                             //
    float sum;                                      // temp accumulation buffer

    // check validity of params
    if(!in || !out || !kernel) return -1;
    if(dataSizeX <= 0 || kernelSizeX <= 0) return -1;

    // find center position of kernel (half of kernel size)
    kCenterX = (int)kernelSizeX / 2;
    kCenterY = (int)kernelSizeY / 2;

    // init working  pointers
    inPtr = inPtr2 = &in[dataSizeX * kCenterY + kCenterX];  // note that  it is shifted (kCenterX, kCenterY),
    outPtr = out;
    kPtr = kernel;

    // start convolution
    for(i= 0; i < dataSizeY; ++i)                   // number of rows
    {
        // compute the range of convolution, the current row of kernel should be between these
        rowMax = i + kCenterY;
        rowMin = i - dataSizeY + kCenterY;

        for(j = 0; j < dataSizeX; ++j)              // number of columns
        {
            // compute the range of convolution, the current column of kernel should be between these
            colMax = j + kCenterX;
            colMin = j - dataSizeX + kCenterX;

            sum = 0;                                // set to 0 before accumulate

            // flip the kernel and traverse all the kernel values
            // multiply each kernel value with underlying input data
            for(m = 0; m < kernelSizeY; ++m)        // kernel rows
            {
                // check if the index is out of bound of input array
                if(m <= rowMax && m > rowMin)
                {
                    for(n = 0; n < kernelSizeX; ++n)
                    {
                        // check the boundary of array
                        if(n <= colMax && n > colMin)
                            sum += *(inPtr - n) * *kPtr;

                        ++kPtr;                     // next kernel
                    }
                }
                else
                    kPtr += kernelSizeX;            // out of bound, move to next row of kernel

                inPtr -= dataSizeX;                 // move input data 1 raw up
            }

            // convert integer number
            if(sum >= 0) *outPtr = (int)(sum + 0.5f);
//            else *outPtr = (int)(sum - 0.5f)*(-1);
                // For using with image editors like GIMP or others...
            else *outPtr = (int)(sum - 0.5f);
            // For using with a text editor that read ppm images like libreoffice or others...
//            else *outPtr = 0;

            kPtr = kernel;                          // reset kernel to (0,0)
            inPtr = ++inPtr2;                       // next input
            ++outPtr;                               // next output
        }
    }

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION
//////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    if(argc != 5)
    {
        printf("Usage: %s <image-file> <kernel-file> <result-file> <partitions>\n", argv[0]);

        printf("\n\nError, Missing parameters:\n");
        printf("format: ./serialconvolution image_file kernel_file result_file\n");
        printf("- image_file : source image path (*.ppm)\n");
        printf("- kernel_file: kernel path (text file with 1D kernel matrix)\n");
        printf("- result_file: result image path (*.ppm)\n");
        printf("- partitions : Image partitions\n\n");
        return -1;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // READING IMAGE HEADERS, KERNEL Matrix, DUPLICATE IMAGE DATA, OPEN RESULTING IMAGE FILE
    //////////////////////////////////////////////////////////////////////////////////////////////////
    int  partitions, partsize, chunksize, halo, halosize;
    int numtasks, rank;
    int master_partsize=0;
    int offset=0;
    int master_offset=0;
    int position=0;
    long img_position=0;
    double start, tstart=0, tread=0, tcopy=0, tconv=0, tstore=0, treadk=0;
    struct timeval tim;
    //char *buffer;
    int struct_data_r_g_b_size=0;
    int struct_data_size=0;
    FILE *fpsrc=NULL,*fpsrcaux=NULL,*fpdst=NULL;
    ImagenData source=NULL, output=NULL, master_output=NULL;
    // Store number of partitions
    partitions = atoi(argv[4]);

    ////////////////////////////////////////
    //Reading kernel matrix
    ////////////////////////////////////////
    gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
    tstart = start;
    kernelData kern=NULL;
    if ( (kern = leerKernel(argv[2]))==NULL) {
        return -1;
    }
    //The matrix kernel define the halo size to use with the image. The halo is zero when the image is not partitioned.
    if (partitions==1) halo=0;
    else halo = (kern->kernelY/2)*2;
    gettimeofday(&tim, NULL);
    treadk = treadk + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
    printf("0, READKERNEL,%.6lf\n", treadk);

    ////////////////////////////////////////
    //End Reading kernel matrix
    ////////////////////////////////////////

    ////////////////////////////////////////
    //Reading Image Header. Image properties: Magical number, comment, size and color resolution.
    ///////////////////////////////////////
    //Memory allocation based on number of partitions and halo size.
    if ( (source = initimage(argv[1], &fpsrcaux, partitions, halo)) == NULL) {
        return -1;
    }
    fclose(fpsrcaux);
    ////////////////////////////////////////
    //End Reading Image Header. Image properties: Magical number, comment, size and color resolution.
    ///////////////////////////////////////
    struct_data_r_g_b_size = (source->ancho*source->altura / partitions);
    struct_data_size = struct_data_r_g_b_size * sizeof(int) * 3 + sizeof(int)*7 + sizeof(float) + sizeof(char)*300;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);   // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // get current process id

    if(rank == 0){
        double conv_start;
        gettimeofday(&tim, NULL);
        conv_start =  MPI_Wtime();
        char *buffer = malloc(struct_data_size);
        ////////////////////////////////////////
        //Reading Image
        ///////////////////////////////////////
        freeImagestructure(&source);
        gettimeofday(&tim, NULL);
        start = MPI_Wtime();
        //Memory allocation based on number of partitions and halo size.
        if ( (source = initimage(argv[1], &fpsrc, partitions, halo)) == NULL) {
            return -1;
        }
        gettimeofday(&tim, NULL);
        tread =  (MPI_Wtime() - start);
        printf("%d, READFILE,%.6lf\n", rank, tread);
        ////////////////////////////////////////
        //End Reading Image
        ///////////////////////////////////////

        ////////////////////////////////////////
        //Init master output
        ////////////////////////////////////////
        if ( (master_output = duplicateImageData(source, partitions, halo)) == NULL) {
            return -1;
        }
        if ( (output = duplicateImageData(source, partitions, halo)) == NULL) {
            return -1;
        }
        ////////////////////////////////////////
        //End Init master output
        ////////////////////////////////////////

        ////////////////////////////////////////
        //Initialize Image Storing file. Open the file and store the image header.
        ////////////////////////////////////////
        if (initfilestore(master_output, &fpdst, argv[3], &img_position)!=0) {
            perror("Error: ");
            return -1;
        }
        ////////////////////////////////////////
        //End Initialize Image Storing file. Open the file and store the image header.
        ////////////////////////////////////////


        //////////////////////////////////////////////////////////////////////////////////////////////////
        // CHUNK READING
        //////////////////////////////////////////////////////////////////////////////////////////////////
        int c=0;
        partsize  = (source->altura*source->ancho)/partitions;
        while (c < partitions) {
            ////////////////////////////////////////////////////////////////////////////////
            //Reading Next chunk.
            if (c == 0) {
                halosize = halo / 2;
                chunksize = partsize + (source->ancho * halosize);
                offset = 0;
            } else if (c < partitions - 1) {
                halosize = halo;
                chunksize = partsize + (source->ancho * halosize);
                offset = (source->ancho * halo / 2);
            } else {
                halosize = halo / 2;
                chunksize = partsize + (source->ancho * halosize);
                offset = (source->ancho * halo / 2);
            }
            if (readImage(source, &fpsrc, chunksize, halo / 2, &img_position)) {
                return -1;
            }
            if(c < partitions - 1){
                position = 0;
                MPI_Pack(&source->altura, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(&source->ancho, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(source->comentario, 300, MPI_CHAR, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(&source->maxcolor, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(&source->P, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(source->R, struct_data_r_g_b_size, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(source->G, struct_data_r_g_b_size, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(source->B, struct_data_r_g_b_size, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(&partsize, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(&offset, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                MPI_Pack(&halosize, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
                // Position has been incremented by sizeof(float) bytes
                MPI_Send((void*)buffer, struct_data_size, MPI_PACKED, c+1,0, MPI_COMM_WORLD);
            }
            c++;
        }

        ////////////////////////////////////////
        //Init output
        ////////////////////////////////////////
        if ( (master_output = duplicateImageData(source, partitions, halo)) == NULL) {
            return -1;
        }
        ////////////////////////////////////////
        //End Init output
        ////////////////////////////////////////
        if ( duplicateImageChunk(source, master_output, chunksize) ) {
            return -1;
        }
         master_partsize = partsize;
         master_offset = offset;
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // CHUNK CONVOLUTION
        //////////////////////////////////////////////////////////////////////////////////////////////////
        convolve2D(source->R, master_output->R, source->ancho, (source->altura / partitions) + halosize, kern->vkern,
                   kern->kernelX, kern->kernelY);
        convolve2D(source->G, master_output->G, source->ancho, (source->altura / partitions) + halosize, kern->vkern,
                   kern->kernelX, kern->kernelY);
        convolve2D(source->B, master_output->B, source->ancho, (source->altura / partitions) + halosize, kern->vkern,
                   kern->kernelX, kern->kernelY);
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // END CONVOLUTION
        //////////////////////////////////////////////////////////////////////////////////////////////////

        c=0;
        while(c < partitions-1){
            MPI_Recv((void*)buffer,struct_data_size,MPI_PACKED, c+1, MPI_ANY_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            position=0;
            MPI_Unpack(buffer, struct_data_size, &position, &output->altura, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &output->ancho, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, output->comentario, 300, MPI_CHAR, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &output->maxcolor, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &output->P, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, output->R, struct_data_r_g_b_size, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, output->G, struct_data_r_g_b_size, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, output->B, struct_data_r_g_b_size, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &partsize, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &offset, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &halosize, 1, MPI_INT, MPI_COMM_WORLD);


            //////////////////////////////////////////////////////////////////////////////////////////////////
            // CHUNK SAVING
            //////////////////////////////////////////////////////////////////////////////////////////////////
            start =  MPI_Wtime();
            if (savingChunk(output, &fpdst, partsize, offset)) {
                perror("Error: ");
                return -1;
            }
            tstore =  (MPI_Wtime() - start);
            printf("%d, SAVE,%.6lf\n", c+1, tstore);
            //////////////////////////////////////////////////////////////////////////////////////////////////
            // END CHUNK SAVING
            //////////////////////////////////////////////////////////////////////////////////////////////////

           c++;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Master CHUNK SAVING
        //////////////////////////////////////////////////////////////////////////////////////////////////
        tstore=0;
        start =  MPI_Wtime();
        if (savingChunk(master_output, &fpdst, master_partsize, master_offset)) {
            perror("Error: ");
            return -1;
        }
        tstore =  (MPI_Wtime() - start);
        printf("%d, SAVE,%.6lf\n", rank, tstore);
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // End Master CHUNK SAVING
        //////////////////////////////////////////////////////////////////////////////////////////////////
        gettimeofday(&tim, NULL);
        tconv = (MPI_Wtime() - conv_start);
        printf("%d, CONV,%.6lf\n", rank, tconv);
        freeImagestructure(&master_output);
        fclose(fpsrc);
        fclose(fpdst);
    }else{
        char *buffer = malloc(struct_data_size);
        if(partitions != 1){
            MPI_Recv((void*)buffer,struct_data_size, MPI_PACKED,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            position=0;

            MPI_Unpack(buffer, struct_data_size, &position, &source->altura, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &source->ancho, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, source->comentario, 300, MPI_CHAR, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &source->maxcolor, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &source->P, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, source->R, struct_data_r_g_b_size, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, source->G, struct_data_r_g_b_size, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, source->B, struct_data_r_g_b_size, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &partsize, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &offset, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, struct_data_size, &position, &halosize, 1, MPI_INT, MPI_COMM_WORLD);


            ////////////////////////////////////////
            //Init output
            ////////////////////////////////////////
            if ( (output = duplicateImageData(source, partitions, halo)) == NULL) {
                return -1;
            }
            ////////////////////////////////////////
            //End Init output
            ////////////////////////////////////////
            chunksize = partsize + (source->ancho * halosize);
            if ( duplicateImageChunk(source, output, chunksize) ) {
                return -1;
            }

            //////////////////////////////////////////////////////////////////////////////////////////////////
            // CHUNK CONVOLUTION
            //////////////////////////////////////////////////////////////////////////////////////////////////
            convolve2D(source->R, output->R, source->ancho, (source->altura / partitions) + halosize, kern->vkern,
                       kern->kernelX, kern->kernelY);
            convolve2D(source->G, output->G, source->ancho, (source->altura / partitions) + halosize, kern->vkern,
                       kern->kernelX, kern->kernelY);
            convolve2D(source->B, output->B, source->ancho, (source->altura / partitions) + halosize, kern->vkern,
                       kern->kernelX, kern->kernelY);
            //////////////////////////////////////////////////////////////////////////////////////////////////
            // END CHUNK CONVOLUTION
            //////////////////////////////////////////////////////////////////////////////////////////////////
            //send to master = source + partsize + offset + halosize
            position = 0;
            MPI_Pack(&output->altura, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(&output->ancho, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(output->comentario, 300, MPI_CHAR, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(&output->maxcolor, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(&output->P, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(output->R, struct_data_r_g_b_size, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(output->G, struct_data_r_g_b_size, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(output->B, struct_data_r_g_b_size, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(&partsize, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(&offset, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            MPI_Pack(&halosize, 1, MPI_INT, buffer, struct_data_size, &position, MPI_COMM_WORLD);
            // Position has been incremented by sizeof(float) bytes
            MPI_Send((void*)buffer, struct_data_size, MPI_PACKED, 0,0, MPI_COMM_WORLD);
        }
    }
    freeImagestructure(&source);
    freeImagestructure(&output);
    MPI_Finalize();
    return 0;
}
