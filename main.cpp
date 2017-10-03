// cookbook CCfits demonstration program
//	Astrophysics Science Division,
//	NASA/ Goddard Space Flight Center
//	HEASARC
//	http://heasarc.gsfc.nasa.gov
//	e-mail: ccfits@legacy.gsfc.nasa.gov
//
//	Original author: Ben Dorman


// The CCfits headers are expected to be installed in a subdirectory of
// the include path.

// The <CCfits> header file contains all that is necessary to use both the CCfits
// library and the cfitsio library (for example, it includes fitsio.h) thus making
// all of cfitsio's symbolic names available.

#ifdef _MSC_VER
#include "MSconfig.h" // for truncation warning
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// this includes 12 of the CCfits headers and will support all CCfits operations.
// the installed location of the library headers is $(ROOT)/include/CCfits

// to use the library either add -I$(ROOT)/include/CCfits or #include <CCfits/CCfits>
// in the compilation target.


#include <vector>
#include <opencv2/opencv.hpp>
#include <string>
#include "core/houghUtilities.h"

using namespace cv;
#define DEBUG
#define N 9 //numbers of image, it defines the number of mask
//pinta2 functions is for debug purpose
void pinta2(ImageValChar& val,int Dy,int Dx, int indice);

//void write_im(ImageValChar& val,int Dy,int Dx, int indice);


int main(int argc,char *argv[]){
	int dimX;					//1024 //dimension imagenes HRT de prueba
	int dimY;				//1024
	string nombreImagen;
	//four sets of 9 each displacement image
	//char imageName[] = "./im/im0X.fits"; //   8  images set with a displacement below 15% of  solar disc radius
	//char imageName[] = "./imF03/im0X.fits";//11 images set with a displacement around 0.3% of  solar disc radius
	//char imageName[] = "./imF1/im0X.fits";//10  images set with a displacement up to 20% of  solar disc radius
	//char imageName[] = "./imF2/im0X.fits";//  10  images set with a displacement up to 40% of  solar disc radius

	//Test_SolarC
	//char imageName[] = "./Test_SolarC/im0X.fits"; // 9  images 1080x1080 radio 375 desplazamiento maximo 300 entrono al centro, centro 540,540
	vector <ImageValInt> datacube;

   /* float Kernel[3][3] = {
                          {1/9.0, 1/9.0, 1/9.0},
                          {1/9.0, 1/9.0, 1/9.0},
                          {1/9.0, 1/9.0, 1/9.0}
                         };*/

    float Kernel[3][3] = {
                              {1/16.0, 2/16.0, 1/16.0},
                              {2/16.0, 4/16.0, 2/16.0},
                              {1/16.0, 2/16.0, 1/16.0}
                             };

	ImageValChar tmp (dimX*dimY); //to store masks

	float rand_parameter=0.5;//It defines the size of the solar limbs coordinates subset


			int numero;
		    vector<string> dir=read_directory(argv[1],numero);
			//string filename=argv[1]+dir[0];
			//vector<vector<int>> disp;
			//disp=read_disp(filename);//this function gets displacements
			cout << "!!!Creating data cube of images.....   " <<endl<< endl;
	// read images from fits files into vector data cube 3D objects
	for(unsigned int i = 0; i < 1; i++) {

		//imageName[8] = 48 + i;//for "./im/im0X.fits" set
		//imageName[11] = 48 + i;
		//imageName[17] = 48 + i;//for "./Test_SolarC/im0X.fits" set
		//nombreImagen = imageName;
		string tmp = argv[1];
		string tmp2 = dir[i];
		nombreImagen = tmp + "/" + tmp2;
		cout << nombreImagen  << endl;
		datacube.push_back(readImageFit(nombreImagen,dimX,dimY));
		//datacube.push_back(readImageFit("./im/im00.fits",dimX,dimY));


		//Calculate image gradient

		ImageValInt ima=datacube[i];


		//ima2=rot90(ima,dimX,dimY);
		//valarray<double> ImFil=median_filter(ima,Kernel);
		//ImageValDouble valorG = gradient(ImFil);
		ImageValLong valorG = gradient(ima,dimX,dimY);
		cout<< "maximo del gradiente    "<< valorG.max() <<endl<<endl;
		ImageValChar valorG8 = escalado8(valorG);

		//pinta2(valorG8,dimX,dimY,i+1);
		int umbral = otsu_th(valorG8,dimX,dimY);
		cout<< "umbral    "<< umbral <<endl<<endl;
		binarizar(valorG8,umbral);
		ImageValInt ones = findones(valorG8,dimX,dimY); //find solar limbs coordinates
		ImageValInt random_ones = randomizer(ones, rand_parameter);
		//hough function with initial parameters
		int Xc,Yc;
		int Rinicial=963;
		int Range; //range is square side size of center searching
		Xc=dimX/2;
		Yc=dimY/2;
		Range=(int)Xc*30/100;

		//coarse pass (pixel precision) and find the centers around the CCD centers in ROI (600*600) square

		//ImageValFloat matrix = hough(random_ones, 963.8, 1, 1020.68, 1021.75, 600);
		//test SOlar C
		//coarse pass (pixel precision) and find the centers around the CCD centers in ROI (300*300) square
		//ImageValFloat matrix = hough(random_ones, 376, 1, 540, 540, 300);
		ImageValFloat matrix = hough(random_ones, Rinicial, 1, Xc, Yc, Range);

#ifdef DEBUG

		//		cout<< "image name"<< nombreImagen<<endl;
		//		cout << "Otsu threhol : " << umbral << endl;
		//		cout << "Size Ones: " << ones.size()/2 << endl;
		//		cout << "Size Random_Ones: " << random_ones.size()/2 << endl;
		//		cout<< "randon _parameter"<< rand_parameter<<endl;
		//		write_im(valorG8,dimX,dimY, i);
		//		ImageValChar valor = escalado8(ima);
		//		pinta2(valor,dimX,dimY,i);
		//		waitKey(0);
		//		valor = escalado8(ima2);
		//		pinta2(valor,dimX,dimY,i+1);
		//find max in matrix vector

		int largest = matrix[0];
		int indice;
		for(int i=0; i < matrix.size(); i+=4){
			if(i%4==0){
				//cout << endl;
			}
			if (largest<matrix[i]) {   // or: if (comp(*largest,*first)) for version (2)
				largest=matrix[i];
				indice=i;
			}
			//cout << matrix[i] << "    " << matrix[i+1] << "    "  << matrix[i+2] << "    "  << matrix[i+3] << endl;
		}
		//waitKey(0);
		cout << matrix[indice] << "    " << matrix[indice+1] << "    "  << matrix[indice+2] << "    "  << matrix[indice+3] << endl;

	}
	cout<< "finnn"<<endl;

#endif

	return 0;
}

void pinta2(ImageValChar& val,int Dy,int Dx, int indice){

	Mat im(Dy, Dx, CV_8U, Scalar(0));  //Es un tipo de dato de 4 bytes 32S


	//Se pone primero el eje Y y despues el eje XCV_64F
	for (int y=0; y<Dy; y++){

		for (int x=0; x<Dx; x++){
			//cout << " y   x  : "  << y*Dx + x << "  " << x << endl;
			im.at<uchar>(y,x) = val[y*Dx + x];
		}
	}
	char imageName[] = "imX.jpg";
	imageName[2] = 48 + indice;
	imwrite(imageName, im);

	namedWindow(imageName, CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
	imshow(imageName, im);
}

