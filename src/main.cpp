/*
*   Example on how to run the program
*   ./main --directory ~/reconstruction/Voxel-Carving/assets --voxdim 128
*   images have just names 00.jpg->99.jpg or 000.jpg->999.jpg
*   images must be in this format to be lined up correctly at run time
*   Evan O'Keeffe
*   10324289
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>


#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkPLYWriter.h>
#include <vtkFloatArray.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkMarchingCubes.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

using namespace std;
using namespace boost::program_options;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

int IMG_WIDTH = 1280;
int IMG_HEIGHT = 960;
int VOXEL_DIM = 128;
int VOXEL_SIZE = VOXEL_DIM*VOXEL_DIM*VOXEL_DIM;
int VOXEL_SLICE = VOXEL_DIM*VOXEL_DIM;
const int OUTSIDE = 0;
string directory = "";
string output = "";
bool render_it = false;

struct voxel {
    float xpos;
    float ypos;
    float zpos;
    float res;
    float value;
};

struct coord {
    int x;
    int y;
};

struct startParams {
    float startX;
    float startY;
    float startZ;
    float voxelWidth;
    float voxelHeight;
    float voxelDepth;
};

struct camera {
    cv::Mat Image;
    cv::Mat P;
    cv::Mat K;
    cv::Mat R;
    cv::Mat t;
    cv::Mat Silhouette;
};

void exportModel(std::string filename, vtkSmartPointer<vtkCleanPolyData> polyData) {
    
    /* exports 3d model in ply format */
    vtkSmartPointer<vtkPLYWriter> plyExporter = vtkSmartPointer<vtkPLYWriter>::New();
    plyExporter->SetFileName(filename.c_str());
    plyExporter->SetInputConnection(polyData->GetOutputPort());
    plyExporter->Update();
    plyExporter->Write();
}

coord project(camera cam, voxel v) {
    
    coord im;
    
    /* project voxel into camera image coords */
    float z =   cam.P.at<float>(2, 0) * v.xpos +
                cam.P.at<float>(2, 1) * v.ypos +
                cam.P.at<float>(2, 2) * v.zpos +
                cam.P.at<float>(2, 3);
    
    im.y =    (cam.P.at<float>(1, 0) * v.xpos +
               cam.P.at<float>(1, 1) * v.ypos +
               cam.P.at<float>(1, 2) * v.zpos +
               cam.P.at<float>(1, 3)) / z;
    
    im.x =    (cam.P.at<float>(0, 0) * v.xpos +
               cam.P.at<float>(0, 1) * v.ypos +
               cam.P.at<float>(0, 2) * v.zpos +
               cam.P.at<float>(0, 3)) / z;
    
    return im;
}

void renderModel(float fArray[], startParams params,std::string filename) {
    
    /* create vtk visualization pipeline from voxel grid (float array) */
    vtkSmartPointer<vtkStructuredPoints> sPoints = vtkSmartPointer<vtkStructuredPoints>::New();
    sPoints->SetDimensions(VOXEL_DIM, VOXEL_DIM, VOXEL_DIM);
    sPoints->SetSpacing(params.voxelDepth, params.voxelHeight, params.voxelWidth);
    sPoints->SetOrigin(params.startZ, params.startY, params.startX);
    sPoints->SetScalarTypeToFloat();
    
    vtkSmartPointer<vtkFloatArray> vtkFArray = vtkSmartPointer<vtkFloatArray>::New();
    vtkFArray->SetNumberOfValues(VOXEL_SIZE);
    vtkFArray->SetArray(fArray, VOXEL_SIZE, 1);
    
    sPoints->GetPointData()->SetScalars(vtkFArray);
    sPoints->Update();
    
    /* create iso surface with marching cubes algorithm */
    vtkSmartPointer<vtkMarchingCubes> mcSource = vtkSmartPointer<vtkMarchingCubes>::New();
    mcSource->SetInputConnection(sPoints->GetProducerPort());
    mcSource->SetNumberOfContours(1);
    mcSource->SetValue(0,0.5);
    mcSource->Update();
    
    /* recreate mesh topology and merge vertices */
    vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanPolyData->SetInputConnection(mcSource->GetOutputPort());
    cleanPolyData->Update();
    
    /* usual render stuff */
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(.45, .45, .9);
    renderer->SetBackground2(.0, .0, .0);
    renderer->GradientBackgroundOn();
    
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);
    
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(cleanPolyData->GetOutputPort());
    
    exportModel(filename,cleanPolyData);

    if(render_it)
    {
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        
        /* visible light properties */
        actor->GetProperty()->SetSpecular(0.15);
        actor->GetProperty()->SetInterpolationToPhong();
        renderer->AddActor(actor);

        renderWindow->Render();
        interactor->Start();
    }
}

void carve(float fArray[], startParams params, camera cam) {
    
    cv::Mat silhouette, distImage;
    cv::Canny(cam.Silhouette, silhouette, 0, 255);
    cv::bitwise_not(silhouette, silhouette);
    cv::distanceTransform(silhouette, distImage, CV_DIST_L2, 3);
    
    for (int i=0; i<VOXEL_DIM; i++) {
        for (int j=0; j<VOXEL_DIM; j++) {
            for (int k=0; k<VOXEL_DIM; k++) {
                
                /* calc voxel position inside camera view frustum */
                voxel v;
                v.xpos = params.startX + i * params.voxelWidth;
                v.ypos = params.startY + j * params.voxelHeight;
                v.zpos = params.startZ + k * params.voxelDepth;
                v.value = 1.0f;
                
                coord im = project(cam, v);
                float dist = -1.0f;
                                
                /* test if projected voxel is within image coords */
                if (im.x > 0 && im.y > 0 && im.x < IMG_WIDTH && im.y < IMG_HEIGHT) {
                    dist = distImage.at<float>(im.y, im.x);
                    if (cam.Silhouette.at<uchar>(im.y, im.x) == OUTSIDE) {
                        dist *= -1.0f;
                    }
                }
                
                if (dist < fArray[i*VOXEL_SLICE+j*VOXEL_DIM+k]) {
                    fArray[i*VOXEL_SLICE+j*VOXEL_DIM+k] = dist;
                }
                
            }
        }
    }
    
}

int main(int argc, char* argv[]) {
    
    /* acquire camera images, silhouettes and camera matrix */
    std::vector<camera> cameras;

    po::options_description desc("Showing Options");
    desc.add_options()
        ("help,h", "produce help message")
        ("directory", po::value<string>(&directory), "directry to read")
        ("output",po::value<string>(&output), "output file name")
        ("render",po::value<bool>(&render_it), "render the file with vtk,default is false")
        ("voxdim", po::value<int>(&VOXEL_DIM), "set voxel dimensions")
        ;

    po::variables_map vm;
    // parse regular options
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm); 

    if (fs::is_directory(directory)) 
    {
        cout << "Reading Images from:" << directory << endl;
    } 
    else 
    {
        cout << "No directory given" << endl;
        cout << desc << endl;
        return -1;
    }

    if (vm.count("voxdim")) 
    {
        cout << "Using Voxel Dimension:" << VOXEL_DIM << endl;
    } 
    else 
    {
        cout << desc << endl;
        return -1;
    }

    VOXEL_SIZE = VOXEL_DIM*VOXEL_DIM*VOXEL_DIM;
    VOXEL_SLICE = VOXEL_DIM*VOXEL_DIM;

    std::cout << "Looking for:" << directory << std::endl;
    std::cout << "Looking for:" << directory+"viff.xml" << std::endl;

    cv::FileStorage fs(directory+"viff.xml", cv::FileStorage::READ);
    
    if(!fs.isOpened())
    {
        std::cerr << "No viff.xml file or no permission to access in :" << directory << std::endl;
        return 1;
    }
    else
    {
        std::cout << "viff.xml file found and open" << std::endl;
    }

    fs::directory_iterator it(directory), eod;

    std::vector<string> files;
    BOOST_FOREACH(fs::path const &file_path, std::make_pair(it, eod))
    {
        files.push_back(file_path.string());
    }

    std::sort(files.begin(), files.end());

    int i=0;
    BOOST_FOREACH(std::string const &simg, files)
    {
        //std::string simg = file_path.string();

        if(simg.find(".jpg") == std::string::npos
            &&
           simg.find(".JPG") == std::string::npos
        )
        {
            // essentially if .jpg or .jpg
            continue;
        }

        /* camera image */
        cv::Mat img = cv::imread(simg);
        
        /* silhouette */
        cv::Mat silhouette;
        cv::cvtColor(img, silhouette, CV_BGR2HSV);
        cv::inRange(silhouette, cv::Scalar(0, 0, 30), cv::Scalar(255,255,255), silhouette);
        
        /* camera matrix */
        std::stringstream smat;
        smat << "viff" << std::setfill('0') << std::setw(3) << i++ << "_matrix";
        cv::Mat P;
        fs[smat.str()] >> P;

        /* decompose proj matrix to camera and rotation matrix and translation vectors */
        cv::Mat K, R, t;

        std::cout << smat.str() << std::endl;
        std::cout << P << std::endl;

        cv::decomposeProjectionMatrix(P, K, R, t);
        K = cv::Mat::eye(3, 3, CV_32FC1);
        K.at<float>(0,0) = 1680.2631413061415; /* fx */
        K.at<float>(1,1) = 1676.1202869984309; /* fy */
        K.at<float>(0,2) = 621.59194200994375; /* cx */
        K.at<float>(1,2) = 467.7223229477861; /* cy */
        
        camera c;
        c.Image = img;
        c.P = P;
        c.K = K;
        c.R = R;
        c.t = t;
        c.Silhouette = silhouette;
                
        cameras.push_back(c);
    }
    
    /* bounding box dimensions of squirrel */
    float xmin = -6.21639, ymin = -10.2796, zmin = -14.0349;
    float xmax = 7.62138, ymax = 12.1731, zmax = 12.5358;
            
    float bbwidth = std::abs(xmax-xmin)*1.15;
    float bbheight = std::abs(ymax-ymin)*1.15;
    float bbdepth = std::abs(zmax-zmin)*1.05;
    
    startParams params;
    params.startX = xmin-std::abs(xmax-xmin)*0.15;
    params.startY = ymin-std::abs(ymax-ymin)*0.15;
    params.startZ = 0.0f;
    params.voxelWidth = bbwidth/VOXEL_DIM;
    params.voxelHeight = bbheight/VOXEL_DIM;
    params.voxelDepth = bbdepth/VOXEL_DIM;
    
    /* 3 dimensional voxel grid */
    float *fArray = new float[VOXEL_SIZE];
    std::fill_n(fArray, VOXEL_SIZE, 1000.0f);
    
    std::cout << "carving model now..." << std::endl;
    /* carving model for every given camera image */
    for (int i=0; i<cameras.size(); i++) {
        carve(fArray, params, cameras.at(i));
    }

    /* show example of segmented image */
    cv::Mat original, segmented;
    //cv::resize(cameras.at(1).Image, original, cv::Size(640, 480));
    //cv::resize(cameras.at(1).Silhouette, segmented, cv::Size(640, 480));
    //cv::imshow("Squirrel" , original);
    //cv::imshow("Squirrel Silhouette", segmented);

    std::cout << "Rendering & Saving stl model" << std::endl;
    std::cout << "saving model to output.stl" << std::endl;

    renderModel(fArray, params, output+".stl");
    
    return 0;
}
