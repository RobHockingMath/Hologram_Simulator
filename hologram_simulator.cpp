#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // for std::min, std::max
#include <cstdio>    // for snprintf
#include <omp.h>

#include "minimal_image_io.h" // Provides read_image_as_3d_double and write_image_as_3d_double
#include "raytrace_frame.h"

#define M_PI 3.14159265358979323846




int main() {
    double DY = 40;
    double DX = 40;
	//DY=42.56;
	//DX=48;

    int N_A = 160; // number of hogels horizontally
    int N_B = 160; // number of hogels vertically
    //int N_I = 1980; // angular resolution of each hogel.  
	int N_I = 1920;

	//N_A = 160;
	//N_B = 160;

    int i_r_mid = N_I/2;

    //std::string hogel_folder_name = "HogelImagesHR6_HPara_X/Hogel_Column_";
	//std::string hogel_folder_name = "HogelViews_120HR_Buffer_Brighter/Hogel_Column_";
	//std::string hogel_folder_name = "Oppenheimer_hogels4/Hogel_Column_";
    std::string hogel_folder_name = "Succub_HogelsHR/Hogel_Column_";

    // Load hogels into a vector
    std::vector<unsigned char***> hogels(N_A,nullptr);
    for (int k=0; k<N_A; k++) {
        std::string fname = hogel_folder_name + std::to_string(k) + ".png";
        int w,h;
        if (!read_image_as_3d_uchar(fname.c_str(), w, h, &hogels[k])) {
            std::cerr << "Could not read hogel: " << fname << std::endl;
            return 1;
        }
        // We assume (w,h) == (N_I,N_B). You can add checks:
        // if (w!=N_I || h!=N_B) { std::cerr<<"Hogel size mismatch\n"; return 1; }
    }

    double dx = DX/N_A;
    double dy = DY/N_B;

    int M = 601;

    int J=int(floor(1080*1.5));
    int I=int(floor((1920-1080)*1.5)); // I=840

    // Load T
    unsigned char ***T_img=nullptr;
    int T_w, T_h;
    if(!read_image_as_3d_uchar("sky2_4096_2048.png", T_w, T_h, &T_img)) {
        std::cerr << "Could not read T image.\n";
        return 1;
    }

    double dtheta = 2*M_PI/T_w;
    double dphi = M_PI/T_h;

    std::vector<double> look={0,DX/2,DY/2};
    std::vector<double> up={0,0,1};

    auto norm_vec =[](const std::vector<double>&a){
        double n=0; for (auto v:a)n+=v*v; n=std::sqrt(n);
        std::vector<double>r(a.size());for (size_t i=0;i<a.size();i++) r[i]=a[i]/n;
        return r;
    };
    auto crossp=[](const std::vector<double>&a,const std::vector<double>&b){
        return std::vector<double>{a[1]*b[2]-a[2]*b[1],
                                   a[2]*b[0]-a[0]*b[2],
                                   a[0]*b[1]-a[1]*b[0]};
    };

    double rot_angle = (M_PI/180)*24;

    double EE1[3]={cos(rot_angle),0,sin(rot_angle)};
	double EE3[3]={-sin(rot_angle),0,cos(rot_angle)};
	//double EE2[3]={0,1,0};
	double EE2[3]={EE1[1]*EE3[2]-EE1[2]*EE3[1],EE1[2]*EE3[0]-EE1[0]*EE3[2],EE1[0]*EE3[1]-EE1[1]*EE3[0]};

    up[0]=EE3[0];up[1]=EE3[1];up[2]=EE3[2];

    int N_rays=3;
	#pragma omp parallel for
    for (int l=0; l<M; l++) {
        std::cout << l << std::endl;
        double theta = M_PI/2-2*M_PI/6+4*M_PI/6*l/(M-1);
        //double R=500;
		double R=100;
        std::vector<double> p={R*std::sin(theta),DX/2 - R*std::cos(theta),DY/2};

        //std::vector<double> p={R*EE1[0]*sin(theta)+R*EE2[0]*cos(theta),DX/2 -R*EE1[1]*sin(theta)-R*EE2[1]*cos(theta),DY/2+R*EE1[2]*sin(theta)+R*EE2[2]*cos(theta)};

        //p[0]=front+r_geola*EE1[0]*sin(theta_chimera)+r_geola*EE2[0]*cos(theta_chimera);
		//	p[1]=0+r_geola*EE1[1]*sin(theta_chimera)+r_geola*EE2[1]*cos(theta_chimera);
		//	p[2]=-LL+r_geola*EE1[2]*sin(theta_chimera)+r_geola*EE2[2]*cos(theta_chimera);

        std::vector<double> lp={look[0]-p[0],look[1]-p[1],look[2]-p[2]};
        std::vector<double> e1=norm_vec(lp);
        std::vector<double> e2=norm_vec(crossp(up,e1));
        std::vector<double> e3=crossp(e1,e2);



        double ***U=new double**[I];
        for (int ii=0;ii<I;ii++){
            U[ii]=new double*[J];
            for (int jj=0;jj<J;jj++){
                U[ii][jj]=new double[3];
				U[ii][jj][0]=0;
				U[ii][jj][1]=0;
				U[ii][jj][2]=0;
            }
        }

        for (int ii=0; ii<I; ii++) {
            for (int jj=0; jj<J; jj++) {
                for (int iii=0; iii<N_rays; iii++) {
                    for (int jjj=0; jjj<N_rays; jjj++) {
                        double dii=(double)iii/(N_rays);
                        double djj=(double)jjj/(N_rays);
                        double x_part=(-0.75*DX+1.5*DX*(ii+dii-1)/(I-1));
                        double y_part=(-0.75*DY+1.5*DY*(jj+djj-1)/(J-1));

                        std::vector<double> v={
                            R*e1[0]+x_part*e2[0]+y_part*e3[0],
                            R*e1[1]+x_part*e2[1]+y_part*e3[1],
                            R*e1[2]+x_part*e2[2]+y_part*e3[2]
                        };
                        v=norm_vec(v);

                        double t_f; std::vector<double>C_f;
                        raytrace_frame(p,v,DX,DY,t_f,C_f);

                        double t=-p[0]/v[0];
                        double x=p[1]+t*v[1];
                        double y=p[2]+t*v[2];


                        if (x>=0 && x<=DX && y>=0 && y<=DY && (t_f==-1 || (t_f>0 && t<t_f))) {
                            int i_h = std::min(std::max((int)std::lround(x/dx),0),N_A-1);
                            int j_h = std::min(std::max((int)std::lround(y/dy),0),N_B-1);

                            double t2=(1-p[0])/v[0];
                            double x2=p[1]+t2*v[1];
                            double y2=p[2]+t2*v[2];
                            double r=(x - x2);
                            double r_m = std::tan(0.5*120*M_PI/180.0);
                            double dr=2*r_m/N_I;
                            int i_r=std::min(std::max((int)std::lround((r+r_m)/dr),0),N_I-1);

                            unsigned char ***H_a=hogels[i_h];
                            // hogel is N_BxN_I: j_h is row index, i_r is column index
                            // flip horizontally means col -> (N_I - i_r)
                            int flipped_col=(N_I -1 - i_r);


                            double Rv=double(H_a[flipped_col][j_h][0]);
                            double Gv=double(H_a[flipped_col][j_h][1]);
                            double Bv=double(H_a[flipped_col][j_h][2]);

                            U[ii][jj][0]+=Rv;
                            U[ii][jj][1]+=Gv;
                            U[ii][jj][2]+=Bv;
                        } else {
                            if (t_f>=0) {
                                // C_f in [0..255]? Stub: [0..0]
                                U[ii][jj][0]+=C_f[0];
                                U[ii][jj][1]+=C_f[1];
                                U[ii][jj][2]+=C_f[2];
                            } else {
                                double theta_1 = M_PI+std::atan2(v[1],v[0]);
                                double phi_1 = std::acos(v[2]);
                                int i_b=std::min(std::max((int)std::floor(theta_1/dtheta),0),T_w-1);
                                int j_b=std::min(std::max((int)std::floor(phi_1/dphi),0),T_h-1);

                                double Rv=double(T_img[i_b][j_b][0]);
                                double Gv=double(T_img[i_b][j_b][1]);
                                double Bv=double(T_img[i_b][j_b][2]);

                                U[ii][jj][0]+=Rv;
                                U[ii][jj][1]+=Gv;
                                U[ii][jj][2]+=Bv;
                            }
                        }
                    }
                }
            }
        }
		//normalize
		for (int ii=0; ii<I; ii++) {
            for (int jj=0; jj<J; jj++) {
				U[ii][jj][0]/=(255*N_rays*N_rays);
				U[ii][jj][1]/=(255*N_rays*N_rays);
				U[ii][jj][2]/=(255*N_rays*N_rays);
			}
		}

        char outname[256];
        snprintf(outname,256,"Succub_Simulation_HR/%d_space.png",M-1-l);
        if(!write_image_as_3d_double(outname,I,J,U)){
            std::cerr<<"Failed to write "<<outname<<"\n";
        }

		for(int i=0;i<I;i++){
			for(int j=0;j<J;j++){
				delete[] U[i][j];
			}
			delete[] U[i];
		}
		delete[] U;


    }

}