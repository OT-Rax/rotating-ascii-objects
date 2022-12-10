#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#define PI (3.141592653589793)

const int screen_width = 50;
const int screen_height = 50;
const double theta_spacing = 0.07; //theta angle sampling
const double phi_spacing = 0.02; //phi angle sampling
const double R1 = 1; //torus building circles' radius (changes torus size)
const double R2 = 2; //first circle center x coordinate (center = (R2, 0, 0))
const double K2 = 5; //torus depth (distance from screen, summed on z axis)
const double K1 = screen_width*K2*3/(8*(R1+R2));

void renderFrame(double alpha, double beta){
	double cosAlpha = cos(alpha), sinAlpha = sin(alpha);
	double cosBeta = cos(beta), sinBeta = sin(beta);

	char output[screen_width][screen_height];
	memset(output, ' ', screen_height*screen_width*sizeof(char));

	double zbuffer[screen_width][screen_height];
	memset(zbuffer, 0.0, screen_height*screen_width*sizeof(double));

	for (double theta=0; theta < 2*PI; theta+=theta_spacing){
		double cosTheta = cos(theta), sinTheta = sin(theta);
		
		for(double phi=0; phi <2*PI; phi+=phi_spacing){
			double cosPhi = cos(phi), sinPhi = sin(phi);

			double circlex = R2 + R1*cosTheta;
			double circley = R1*sinTheta;

			double x = circlex*(cosBeta*cosPhi + sinAlpha*sinBeta*sinPhi) - circley*cosAlpha*sinBeta;
			double y = circlex*(sinBeta*cosPhi - sinAlpha*cosBeta*sinPhi) + circley*cosAlpha*cosBeta;
			double z = K2 + cosAlpha*circlex*sinPhi + circley*sinAlpha;
			double zInverse = 1/z;

			int xp = (int) (screen_width/2 + K1*zInverse*x);
			int yp = (int) (screen_height/2 - K1*zInverse*y);

			double L = cosPhi*cosTheta*sinBeta - cosAlpha*cosTheta*sinPhi - sinAlpha*sinTheta + cosBeta*(cosAlpha*sinTheta - cosTheta*sinAlpha*sinPhi);

			if(L>0){
				if(zInverse > zbuffer[xp][yp]){
					int index = L*8;
					output[xp][yp] = ".,-~:;=!*#$@"[index];
				}
			}
		}
	}

	//printf("\e[1;1H\e[2J");
	printf("\x1b[H");
	for (int j = 0; j < screen_height; j++){
		for (int i = 0; i <screen_width; i++){
			putchar(output[i][j]);
		}
		putchar('\n');
	}
}


int main(int argc, char* argv[]){
	
	double fps = 60;
	double sleepingTime = 1/fps;
	double alpha = 0;
	double beta = 0;
	double alphaPerSecond= PI/8; 
	double betaPerSecond= PI/4;
	double alphaPerFrame = alphaPerSecond/fps;
	double betaPerFrame = betaPerSecond/fps;

	while(1){
		renderFrame(alpha, beta);	
		alpha += alphaPerFrame;
		beta += betaPerFrame;
		sleep(sleepingTime);
	}
}
