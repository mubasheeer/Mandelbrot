#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include<pthread.h>
#include<sys/time.h>
#include<time.h>

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image( );
void* threading (void* threadargs);

struct bitmap *bm;

struct position{
	int start_point;
	int end_point;
	int height_iterator;
};

struct position my_thread_position[512];
struct position* my_thread_pointer = my_thread_position;

int max;
double xmin;
double xmax;
double ymin;
double ymax;

int no_threads = 2 ;


void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

int main( int argc, char *argv[] )
{
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	max = 1000;

	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:n:o:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'n':
				no_threads = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s\n",xcenter,ycenter,scale,max,outfile);

	// Create a bitmap of the appropriate size.
	bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	xmin=xcenter-scale;
	xmax=xcenter+scale;
	ymin=ycenter-scale;
	ymax=ycenter+scale;
	double time_spent = 0.0;
	clock_t begin = clock();
	compute_image();
	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;

	printf("Time elpased is %f seconds\n", time_spent);

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image()
{
	int i,j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);
	int adjust_width = width/no_threads;
	//printf("adjust width=%d\n",adjust_width );
	pthread_t threads[no_threads];

	// For every pixel in the image...
	my_thread_pointer[0].start_point= 0;
	my_thread_pointer[0].end_point = adjust_width;
	my_thread_pointer[0].height_iterator = 0;

	for(i=1;i<no_threads;i++){
		my_thread_pointer[i].start_point = my_thread_pointer[i-1].end_point;
		my_thread_pointer[i].end_point = my_thread_pointer[i].start_point + adjust_width;
	}

	for(j=0;j<height;j++) {
		for(i=0;i<no_threads;i++){
			my_thread_pointer[i].height_iterator = j;
			pthread_create(&threads[i],NULL,threading,&my_thread_pointer[i]);
		}
		for(i=0;i<no_threads;i++){
			pthread_join(threads[i],NULL);
		}
	}
}
void* threading (void* threadargs){

	struct position* thread_struct;
	thread_struct = (struct position*) threadargs;
	int i;
	int width = bitmap_width(bm);
	int height = bitmap_height(bm);
	int j = thread_struct->height_iterator;


	for(i=thread_struct->start_point;i<thread_struct->end_point;i++) {

		// Determine the point in x,y space for that pixel.
		double x = xmin + i*(xmax-xmin)/width;
		double y = ymin + j*(ymax-ymin)/height;

		// Compute the iterations at that point.
		int iters = iterations_at_point(x,y,max);

		// Set the pixel in the bitmap.
		bitmap_set(bm,i,j,iters);
	}

	return NULL;

}
/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}
