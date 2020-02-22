#ifndef NDT_MCL_3D_HPP_
#define NDT_MCL_3D_HPP_

#include <pcl/conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/io/io.h>
#include <cstdio>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include <fstream>
#include <ndt_map/ndt_map.h>
#include <ndt_map/ndt_cell.h>
#include <ndt_map/pointcloud_utils.h>
#include "ndt_mcl/ParticleFilter3D.h"
#include <tf_conversions/tf_eigen.h>

/**
 * NDT MCL - Class implementation
 */
class NDTMCL3D{
    public:
	lslgeneric::NDTMap map; 		///<This is my map 
	ParticleFilter3D pf; 						///<This is the particle filter
	double resolution;
	double resolution_sensor;
	int counter;
	double zfilt_min;
	bool forceSIR;
        double SIR_varP_threshold;
        int SIR_max_iters_wo_resampling;
        std::vector<double> motion_model, motion_model_offset;

        /**
	 * Constructor
	 */
	NDTMCL3D(bool load_map_from_PCD,std::string mapName,double size_x_meters,double size_y_meters,double size_z_meters,double pose_init_x,double pose_init_y,double map_resolution, lslgeneric::NDTMap &nd_map, double zfilter):
          map(new lslgeneric::LazyGrid(map_resolution)), SIR_varP_threshold(0.006), SIR_max_iters_wo_resampling(25)
	{
            isInit = false;
	    forceSIR = false;
	    resolution=map_resolution;
	    resolution_sensor = resolution;
	    counter = 0;
	    zfilt_min = zfilter;
	    sinceSIR = 0;
		if(load_map_from_PCD){
			map.initialize(pose_init_x,pose_init_y,0,size_x_meters,size_y_meters,size_z_meters);
			pcl::PointCloud<pcl::PointXYZ> pointcloud;
			pcl::io::loadPCDFile (mapName, pointcloud);
			for(unsigned int i=0;i<pointcloud.points.size();i++){
				Eigen::Vector3d m;
				m[0]=pointcloud.points[i].x;
				m[1]=pointcloud.points[i].y;
				m[2]=pointcloud.points[i].z;
				if(m[2]>zfilter){
					Eigen::Matrix3d cov = Eigen::Matrix3d::Identity(3,3);//単位行列　各要素に相関なしとする
					unsigned int nump = 1;//ボクセル内点数1点とする
					//fprintf(stderr,"-%d-",nump);
					map.addDistributionToCell(cov, m,nump);
				}
			}
		}
		else{
		///First, lets make our target map match the given map
	    ///This is done because we want (possibly) lower resolution target map
			double cx,cy,cz;
			if(!nd_map.getCentroid(cx, cy, cz)){
				fprintf(stderr,"Centroid NOT Given-abort!\n");
				exit(1);
			}else{
				fprintf(stderr,"Centroid(%lf,%lf,%lf)\n",cx,cy,cz);
			}

			double wx,wy,wz;
			if(!nd_map.getGridSizeInMeters(wx, wy, wz)){
				fprintf(stderr,"Grid size NOT Given-abort!\n");
				exit(1);
			}else{
				fprintf(stderr,"GridSize(%lf,%lf,%lf)\n",wx,wy,wz);
			}

			map.initialize(cx,cy,cz,wx,wy,wz);

			std::vector<lslgeneric::NDTCell*> ndts;
			ndts = nd_map.getAllCells();
			fprintf(stderr,"NDT MAP with %d components",(int)ndts.size());
			for(unsigned int i=0;i<ndts.size();i++){
				Eigen::Vector3d m = ndts[i]->getMean();	
				if(m[2]>zfilter){
					Eigen::Matrix3d cov = ndts[i]->getCov();
					unsigned int nump = ndts[i]->getN();
					//fprintf(stderr,"-%d-",nump);
					map.addDistributionToCell(cov, m,nump);
				}
			}

		}

	    
	    
		

            motion_model.push_back(0.05);
            motion_model.push_back(0.05);
            motion_model.push_back(0.02);
            motion_model.push_back(0.01);
            motion_model.push_back(0.01);
            motion_model.push_back(0.02);
    
            motion_model.push_back(0.05);
            motion_model.push_back(0.1);
            motion_model.push_back(0.02);
            motion_model.push_back(0.01);
            motion_model.push_back(0.01);
            motion_model.push_back(0.02);

            motion_model.push_back(0.01);
            motion_model.push_back(0.01);
            motion_model.push_back(0.1);
            motion_model.push_back(0.001);
            motion_model.push_back(0.001);
            motion_model.push_back(0.001);

            motion_model.push_back(0.001);
            motion_model.push_back(0.01);
            motion_model.push_back(0.01);
            motion_model.push_back(0.1);
            motion_model.push_back(0.01);
            motion_model.push_back(0.01);

            motion_model.push_back(0.01);
            motion_model.push_back(0.001);
            motion_model.push_back(0.01);
            motion_model.push_back(0.01);
            motion_model.push_back(0.1);
            motion_model.push_back(0.01);

            motion_model.push_back(0.1);
            motion_model.push_back(0.01);
            motion_model.push_back(0.001);
            motion_model.push_back(0.01);
            motion_model.push_back(0.01);
            motion_model.push_back(0.1);

            motion_model_offset.push_back(0.002);
            motion_model_offset.push_back(0.002);
            motion_model_offset.push_back(0.002);
            motion_model_offset.push_back(0.001);
            motion_model_offset.push_back(0.001);
            motion_model_offset.push_back(0.001);
	}

	/**
	 * Intialize filter to some pose and with some number of particles 
	 *
	 */
	void initializeFilter(double x, double y, double z, double roll, double pitch, double yaw,
		double x_e, double y_e, double z_e, double roll_e, double pitch_e, double yaw_e, 
		unsigned int numParticles)
	{
	    pf.initializeNormalRandom(numParticles, x,y,z,roll,pitch, yaw, x_e,y_e,z_e,roll_e,pitch_e,yaw_e);
	}

    void predict(Eigen::Affine3d Tmotion);

	void updateAndPredict(Eigen::Affine3d Tmotion, pcl::PointCloud<pcl::PointXYZ> &cloud);

	void updateAndPredictEff(Eigen::Affine3d Tmotion, pcl::PointCloud<pcl::PointXYZ> &cloud, double subsample_level);
	Eigen::Vector3d getMeanVector(){
	    return (pf.getMean().translation());
	    //mcl::pose m = pf.getDistributionMean(true);
	    //Eigen::Vector3d mm; mm<<m.x,m.y,m.a;
	    //return mm;
	}
	Eigen::Affine3d getMean(){
	    return (pf.getMean());
	    //mcl::pose m = pf.getDistributionMean(true);
	    //Eigen::Vector3d mm; mm<<m.x,m.y,m.a;
	    //return mm;
	}
  

    private:
	bool isInit;
	int sinceSIR;
	double getDoubleTime()
	{
	    struct timeval time;
	    gettimeofday(&time,NULL);
	    return time.tv_sec + time.tv_usec * 1e-6;
	}

};

#endif

