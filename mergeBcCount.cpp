#include<iostream>
#include <fstream>
#include <string>
#include <stdint.h>
#include <map>
#include <vector>
#include <stdlib.h>
#include <hdf5.h>
#include<time.h>

using namespace::std;
void line_split(string line_info,char sep,vector<string> &elements){
  elements.clear();
  string element;
  for(string::size_type ix=0;ix!=line_info.size();ix++){
    if(line_info[ix]!=sep){
      element+=line_info[ix];
    }else{
      elements.emplace_back(element);
      element="";
    }
  }
  elements.emplace_back(element);
}

void errlog(char *pbuf)
{
    char tbuf[32]={0};
    time_t lt1;
    time(&lt1);
    struct tm * tinfo=localtime(&lt1);
    strftime(tbuf, 32, "[%Y-%m-%d %H:%M:%S]\t", tinfo);
    ofstream errout("errcode.log");
    errout<<tbuf<<pbuf<<'\n';
}

int main(int argc,char* argv[]){
	if(argc<4){
	  cerr<<"Arguments error!"<<endl;
	  cerr<<argv[0]<<"\t<mask file>  <input barcode reads count files(a string separated by comma>  <output>"<<endl;
      char buf[64]={0};
      sprintf(buf, "SAW-A30001\t Arguments num error");
      errlog(buf);
	  exit(1);
	}
	string maskFile=argv[1];
  uint32_t range = 0;
  uint32_t offSetInBin=UINT32_MAX;
	if(maskFile.rfind(".h5")==maskFile.size()-3) {
	  string h5FilePath = argv[1];
	  hid_t mask_h5_id = H5Fopen(h5FilePath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	  if(mask_h5_id==H5I_INVALID_HID) {
		cerr << "Error, fail to open h5 file," << h5FilePath << endl;
        char buf[1024]={0};
        sprintf(buf, "SAW-A30002\tFail to open h5 file %s", argv[1]);
        errlog(buf);
		exit(1);
	  }
	  hid_t dataset_id = H5Dopen2(mask_h5_id, "bpMatrix_1", H5P_DEFAULT);
	  hid_t space_id = H5Dget_space(dataset_id);
//  	hid_t plistID = H5Dget_create_plist(dataset_id);
	  int rank = H5Sget_simple_extent_ndims(space_id);
	  hsize_t dims[rank];
	  herr_t status = H5Sget_simple_extent_dims(space_id, dims, NULL);
	  
	  for (int i = 0; i < 2; i++) {
		if(dims[i] > range) {
		  range = dims[i];
		}
	  }
	}else if(maskFile.rfind(".bin")==maskFile.size()-4){
	  ifstream mapReader(maskFile, ios::in | ios::binary);
	  if (! mapReader.is_open()){
		throw invalid_argument("Could not open the file: " + maskFile);
        char buf[1024]={0};
        sprintf(buf, "SAW-A30002\tFail to open bin file %s", argv[1]);
        errlog(buf);
        exit(1);
	  }
//		boost::archive::binary_iarchive ia(mapReader);
//		ia >> bpmap;
	  uint64_t barcodeInt=0;
//		uint64_t position=0;
	  uint32_t posX,posY;
	  while (!mapReader.eof()) {
		mapReader.read((char*)&barcodeInt, sizeof(barcodeInt));
		mapReader.read((char*)&posX, sizeof(posX));
		mapReader.read((char*)&posY, sizeof(posY));
		if(posX>posY){
		  if(posY<offSetInBin){
			offSetInBin=posY;
		  }
		  if(posX>range){
			range=posX;
		  }
		}else{
		  if(posX<offSetInBin){
			offSetInBin=posX;
		  }
		  if(posY>range){
			range=posY;
		  }
		}
//			mapReader.read((char*)&position, sizeof(position));
//			bpmap[barcodeInt] = toPosition1(position);
	  }
	  range-=offSetInBin;
	  mapReader.close();
	}else{
	  cerr<<"Error, only support .h5/.bin suffix mask file"<<endl;
    char buf[1024]={0};
    sprintf(buf, "SAW-A30003\tOnly support .h5/.bin suffix mask file");
    errlog(buf);
	  exit(1);
	}
	  range+=10000; // magic number
	  cout<<"complete cal range"<<endl;
	vector<string> inFiles;
	string p1=argv[2];
	if(p1.find(",")!=string::npos){
		line_split(p1,',',inFiles);
	}else{
		inFiles.push_back(p1);
	}
	//ifstream list(argv[1]);
	string file;
	
	uint16_t** posCount=new uint16_t*[range];
	for(int i=0;i<range;i++){
		posCount[i]=new uint16_t[range];
		for(int j=0;j<range;j++){
			posCount[i][j]=0;
		}
	}
	ofstream out(argv[3]);
	if(!out){
	  cerr<<"Error, fail to open the file,"<<argv[3]<<endl;
    char buf[1024]={0};
    sprintf(buf, "SAW-A30004\tFail to create the file %s", argv[3]);
    errlog(buf);
	  exit(1);
	}
	map<uint64_t,int> overflowCounts;
	//while(getline(list,file)){
	for(vector<string>::iterator ix=inFiles.begin();ix!=inFiles.end();ix++){
		ifstream inFile((*ix).c_str());
		if(!inFile){
		  cerr<<"Error, fail to open the file,"<<*ix<<endl;
        char buf[1024]={0};
        sprintf(buf, "SAW-A30005\tFail to open the file %s", (*ix).c_str());
        errlog(buf);
		  exit(1);
		}
		string line;
		cout<<file<<endl;
		while(getline(inFile,line)){
			int tabNum=0;
			string x,y,count;
			
			for(int i=0;i<line.size();i++){
				if(line[i]=='\t'){
					tabNum++;
				}else{
					switch(tabNum){
						case 0:x+=line[i];break;
						case 1:y+=line[i];break;
						case 2:count+=line[i];break;
						default:break;
					}
				}
			}
			int ix=atoi(x.c_str());
			int iy=atoi(y.c_str());
			int icount=atoi(count.c_str());
			// if(posCount[ix][iy]!=0){
			// 	cerr<<"Warning: found repeat pos,"<<ix<<"\t"<<iy<<endl;
			// }
			if(offSetInBin!=UINT32_MAX){
			  ix-=offSetInBin;
			  iy-=offSetInBin;
			}
			if(ix>=range || iy>=range){
			  cerr<<"Error, out of index, check the mask file whether matched or whether a .h5 converted from .bin"<<endl;
			  cerr<<ix<<"\t"<<iy<<"\t"<<range<<endl;
            char buf[1024]={0};
            sprintf(buf, "SAW-A30006\tOut of index %d %d %d", ix, iy, range);
            errlog(buf);
			  exit(1);
			}
			int sumCount=posCount[ix][iy]+icount;
			if(sumCount>=UINT16_MAX || posCount[ix][iy]==UINT16_MAX){
				uint64_t k=(uint64_t(ix)<<32) | iy;
				if(overflowCounts.find(k)!=overflowCounts.end()){
					overflowCounts[k]+=icount;
				}else{
					overflowCounts[k]=sumCount;
				}
				
				posCount[ix][iy]=UINT16_MAX;
			}else{
				posCount[ix][iy]+=icount;
			}
		}
	}
	
	string outBuf="";
	int num=0;
	for(int i=0;i<range;i++){
		for(int j=0;j<range;j++){
			if(posCount[i][j]>0 && posCount[i][j]<UINT16_MAX){
				num++;
			  if(offSetInBin!=UINT32_MAX){
				out<<i+offSetInBin<<"\t"<<j+offSetInBin<<"\t"<<posCount[i][j]<<"\n";
			  }else {
				out << i << "\t" << j << "\t" << posCount[i][j] << "\n";
			  }
				//outBuf+=to_string(i)+"\t"+to_string(j)+"\t"+to_string(posCount[i][j])+"\n";
                		//if(num>=1000){
                        	//	out<<outBuf;
                        	//	num=0;
                        	//	outBuf="";
               			//}
			}
		}
	}
	//if(outBuf!=""){
	//	cout<<outBuf;
	//}
	for(map<uint64_t,int>::iterator ix=overflowCounts.begin();ix!=overflowCounts.end();ix++){
		int y=(int)ix->first;
		int x=(ix->first)>>32;
		if(offSetInBin!=UINT32_MAX){
		  out << x+offSetInBin << "\t" << y+offSetInBin << "\t" << ix->second << endl;
		}else {
		  out << x << "\t" << y << "\t" << ix->second << endl;
		}
	}
	out.close();
	return 0;
}
