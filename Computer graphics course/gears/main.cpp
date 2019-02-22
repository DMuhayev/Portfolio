#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <tuple>
#include <memory>

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::tuple;
using std::tie;
using std::make_tuple;
using std::shared_ptr;
using std::make_shared;
using std::get;

#include "io.h"
#include "matrix.h"

#include "MyObject.h"

Image Binarisation(const Image& in)
{
	uint i, j;
	uint r, g, b;
	Image out(in.n_rows, in.n_cols);
	
	for (i=0; i<in.n_rows; ++i) {					
		for (j=0; j<in.n_cols; ++j)
		{
				tie (r, g, b) = in(i, j);
				if (r+g+b>130) {r=255; b=255; g=255;} else {r=0; b=0; g=0;}  
				out(i, j) = make_tuple(r, g, b);
		}
	}
	
	return out;
}


tuple<int, vector<shared_ptr<IObject>>, Image>
repair_mechanism(const Image& in, vector<Image> answers)
{
	uint r, g, b, k = 0;
	uint i, j;
	Image out(in.n_rows, in.n_cols);
	
	// binarisation
	out = Binarisation(in);

	for (k=0; k<3; k++)
	{
		answers[k] = Binarisation(answers[k]);
	}
	k=0;
	
	// finding components
	
	vector<uint> labeleq;
	vector<uint> usedlab;
	labeleq.push_back(0);
	
	for (i=0; i<out.n_cols; ++i) 
	{
		tie(r, g, b) = out(0, i);
		if (r!=0) {++k; labeleq.push_back(k); get<0>(out(0, i))= k;}
	}

	for (i=0; i<out.n_rows; ++i) 
	{
		tie(r, g, b) = out(i, 0);
		if (r!=0) {++k; labeleq.push_back(k); get<0>(out(i, 0)) = k;}
	}
			
	uint r1, r2;	

	for (i=1; i<out.n_cols; ++i) {					
		for (j=1; j<out.n_rows; ++j)
		{
				r1 = get<0>(out(j-1, i));
				r2 = get<0>(out(j, i-1));
				tie(r, g, b) = out(j,i);

				if (r==0) continue;
				if ((r1==0) && (r2==0)) {k++; labeleq.push_back(k); get<0>(out(j, i)) = k; continue;}
				if ((r1!=0) && (r2==0)) {get<0>(out(j, i)) = r1; continue;}
				if ((r1==0) && (r2!=0)) {get<0>(out(j, i)) = r2; continue;}
				if (r1==r2) {get<0>(out(j, i)) = r2; continue;}
				if ((r1!=0) && (r2!=0) && (r1!=r2)) {labeleq[r1] = labeleq[r2]; get<0>(out(j, i)) = labeleq[r1];}

		}
	}		


	uint o;
	
	for (i=1; i<=k; i++)
	{
			o=i;
			if (i==labeleq[i]) {continue;}
			while (labeleq[o]!=labeleq[labeleq[o]]) {labeleq[o]=labeleq[labeleq[o]]; o=labeleq[o]; }
	}
	
	for (i=1; i<=k; i++)
	{
		o=0;
		for (j=0; j<usedlab.size(); j++)
		{
				if (labeleq[i]==usedlab[j]) o=1;
		}
		if (o==0) usedlab.push_back(labeleq[i]);
	}
	
	labeleq[0]=usedlab.size();
	
	uint l;  

	for (i=0; i<out.n_rows; ++i) {					
		for (j=0; j<out.n_cols; ++j)
		{
			if (get<0>(out(i, j)) != 0) {l = get<0>(out(i, j)); get<0>(out(i, j)) = labeleq[l];}
		}
	}		

	for (i=0; i<out.n_rows; ++i) {					
		for (j=0; j<out.n_cols; ++j)
		{
			if (get<0>(out(i,j))!=0)
			{
				get<1>(out(i,j)) = ((get<0>(out(i, j))*17)%255)+1;
				get<2>(out(i,j)) = ((get<0>(out(i, j))*7)%255)+1;
			}
		}
	}

//finding gears and axis
	
	auto object_array = vector<shared_ptr<IObject>>();
	Gear gear;
	Axis axis;
	vector<uint> area, areaan; 
	vector<tuple<uint, uint>> border;
	vector<vector<tuple<uint, uint>>> borders;
	vector<vector<tuple<uint, uint>>> ansborders;
	vector<tuple<uint, uint>> centers;
	uint min;

	//finding borders
	for (k=0; k<usedlab.size(); k++)
	{
		for (i=0; i<out.n_rows; ++i) {					
						for (j=0; j<out.n_cols; ++j)
						{	
							r=get<0>(out(i, j));
							if ((r==usedlab[k]) && ( (get<0>(out(i-1, j))!=r) || (get<0>(out(i, j-1))!=r) || (get<0>(out(i+1, j))!=r) || (get<0>(out(i, j+1))!=r) ))
							{
								border.push_back(make_tuple(i, j));
							}
						}	
				}	
		borders.push_back(border);
		border.clear();
	}


// finding centers by distance transform
	vector<tuple<uint, uint>> real_centres;

	for (k=0; k<usedlab.size(); k++)
	{
		for (i=0; i<out.n_rows; i++)	
			{
				for (j=0; j<out.n_cols; j++)
				{
					if (get<0>(out(i,j))==usedlab[k])
						{
							int x = i , y = j;
							int dist = 10000; 
							
							for (l=0; l<borders[k].size(); ++l) 
								{
									int mi, mj;
									tie (mi, mj) = borders[k][l];
									int dd = sqrt( pow(x-mi, 2) + pow(y-mj, 2) );
									
									if (dd < dist) {dist = dd;}
								
								}
							get<1>(out(i,j)) = dist;	
						}
				}
			}
	}
	
	for (k=0; k<usedlab.size(); k++)
	{
		uint max=0;
		uint mi, mj;
		
		for (i=0; i<out.n_rows; i++)	
			{
				for (j=0; j<out.n_cols; j++)
				{
					if (get<0>(out(i,j))==usedlab[k])
						{
							if ( get<1>(out(i,j)) > max ) {max = get<1>(out(i,j)); mi = i; mj = j;}
						}
				}				
			}
		
		real_centres.push_back(make_tuple(mi, mj));	
	}

	// finding area for all segments
	for (k=0; k<usedlab.size(); k++)
	{
		o=0;
			for (i=0; i<out.n_rows; ++i) {					
				for (j=0; j<out.n_cols; ++j)
				{
					r=get<0>(out(i, j));
					if (r==usedlab[k]) {o++;}
				}	
			}	
			
		area.push_back(o);
	}
		
		
// finding area of answers	
	for (k=0; k<3; k++)
	{
		o=0;
		for (i=0; i<answers[k].n_rows; ++i) {					
			for (j=0; j<answers[k].n_cols; ++j)
				{
					r=get<0>(answers[k](i, j));
					if (r!=0) {o++;}
				}	
		}	
		areaan.push_back(o);
	}
	
// centers of answers
	vector<tuple<uint, uint>> anscentres;
	
	
	for (k=0; k<3; k++)
	{
		uint sumi = 0, sumj = 0;
		
		for (i=0; i<answers[k].n_rows; ++i) {					
				for (j=0; j<answers[k].n_cols; ++j)
				{	
					r=get<0>(answers[k](i, j));
					if (r!=0)
					{
						sumi += i;
						sumj += j;
					}
				}	
			
		}	
		
		sumi /= areaan[k];
		sumj /= areaan[k];
		
		anscentres.push_back(make_tuple(sumi, sumj));
	}
	
// 	finding borders of answers
	
	
	for (k=0; k<3; k++)
	{
		for (i=1; i<answers[k].n_rows-1; ++i) {					
						for (j=1; j<answers[k].n_cols-1; ++j)
						{	
							r=get<0>(answers[k](i, j));
							if ((r!=0) && ( (get<0>(answers[k](i-1, j))==0) || (get<0>(answers[k](i, j-1))==0) || (get<0>(answers[k](i+1, j))==0) || (get<0>(answers[k](i, j+1))==0) ))
							{
								border.push_back(make_tuple(i, j));
							}
						}	
				}	
		ansborders.push_back(border);
		border.clear();
	}
	
	uint placei = 0, placej = 0;

// finding centres
	vector<tuple<uint, uint>> centres;
	
	uint ci, cj;
	
	for (k=0; k<usedlab.size(); k++)
	{
		ci=0;
		cj=0;
		for (i=0; i<out.n_rows; ++i) {					
					for (j=0; j<out.n_cols; ++j)
					{	
						r=get<0>(out(i, j));
						if (r==usedlab[k])
						{
							ci += i;
							cj += j;
						}
					}	
			}
			
		ci /= area[k];
		cj /= area[k];
		
		centres.push_back(make_tuple(ci, cj));
	}	

//Is there broken gear? 
	bool broken = 0;
	
	for (k=0; k<usedlab.size(); k++)
	{
		l=0;
		
		int i1, j1, i2, j2;
		i1 = get<0>(centres[k]);
		i2 = get<0>(real_centres[k]);
		j1 = get<1>(centres[k]);
		j2 = get<1>(real_centres[k]); 
		
		if ( (abs(i1-i2)>1)  || (abs(j1-j2)>1) ) {broken=1;} 
	}
	
// finding gears and axis
	min=area[0];
	uint wrong  = 0;
	
	for (k=1; k<area.size(); k++)
	{
			 if (area[k]<min) {min=area[k];}
	}
	
	for (k=0; k<usedlab.size(); k++)
	{
			if (broken!=1)
			{
				if (area[k]==min)  //axis
				{
					
					get<1>(axis.location) = get<0>(centres[k]);
					get<0>(axis.location) = get<1>(centres[k]);
					
					placei = get<1>(axis.location);
					placej = get<0>(axis.location);
					
					object_array.push_back(make_shared<Axis>(axis));
					continue;
				}
				else
				{
					gear.is_broken = 0;
					get<1>(gear.location) = get<0>(centres[k]);
					get<0>(gear.location) = get<1>(centres[k]);
				}
			}
		//gear	
				
			if (broken!=0)
			{
				int i1, j1, i2, j2;
				i1 = get<0>(centres[k]);
				i2 = get<0>(real_centres[k]);
				j1 = get<1>(centres[k]);
				j2 = get<1>(real_centres[k]); 
				
				if ( (abs(i1-i2)>2)  || (abs(j1-j2)>2) ) 
				{
					gear.is_broken = 1;
					get<1>(gear.location) = get<0>(real_centres[k]);
					get<0>(gear.location) = get<1>(real_centres[k]);
					placei = get<1>(gear.location);
					placej = get<0>(gear.location);
					wrong = usedlab[k];
				} 	
				else
				{
					gear.is_broken = 0;
					get<1>(gear.location) = get<0>(centres[k]);
					get<0>(gear.location) = get<1>(centres[k]);
				}
			}	
			
			gear.min_r=10000;
			gear.max_r=0;
			gear.num_cogs=0;
			double rad;
			int sx, sy, si, sj;
			sx = get<1>(gear.location);
			sy = get<0>(gear.location);
			
		// finding radiuses
			for (i=0; i<borders[k].size(); ++i) {
		
							tie (si, sj) = borders[k][i];
							rad = ceil(sqrt(pow((sx-si), 2) + pow((sy-sj), 2)));
							if (rad>gear.max_r) {gear.max_r=rad;};
							if (rad<gear.min_r) {gear.min_r=rad;};
						
			}			
	
		// finding number of cogs
			int radius;
			int distance;
			uint ii, jj;
			
			radius = (gear.max_r+gear.min_r)/2;
			
			for (l=0; l<borders[k].size(); ++l) 
			{
				tie (si, sj) = borders[k][l];
				ii = si;
				jj = sj;
				distance = sqrt(pow((sx-si), 2) + pow((sy-sj), 2));
				if ((distance == radius) || (distance == radius+1)) 
				{
					o=0;
					
					for (i=si-2; i<ii+3; ++i) {					
						for (j=sj-2; j<jj+3; ++j)
						{
							if (get<2>(out(i, j)) == 255) o=1;
						}
					}
							
					if (o==0)		
					{		
					gear.num_cogs++; 
					get<2>(out(si, sj)) = 255;
					}
				}
			}
			
			gear.num_cogs /= 2;
			
			object_array.push_back(make_shared<Gear>(gear));
	}	

//deleting wrong gear
	
	for (i=0; i<out.n_rows; ++i) {					
		for (j=0; j<out.n_cols; ++j)
		{
			if (get<0>(out(i,j)) == wrong )
			{
				out(i, j) = make_tuple(0, 0, 0);
			}
		}
		
	}
	
// finding maximum raduises of answers
	vector<uint> MAR;
	double rad;
	int sx, sy, si, sj;
	
	for (i=0; i<3; i++)
	{
			MAR.push_back(0);
	}
	
	for (k=0; k<3; k++)
	{
		sx = get<0>(anscentres[k]);
		sy = get<1>(anscentres[k]);
		
		for (i=0; i<ansborders[k].size(); ++i) 
		{
			
				tie (si, sj) = ansborders[k][i];
				rad = ceil(sqrt(pow((sx-si), 2) + pow((sy-sj), 2)));
				if (rad>MAR[k]) {MAR[k]=rad;};
							
		}		
	}
	

// finding right answer
    int result_idx = 0;
    vector<bool> Is_placed;
	
	for (k=0; k<3; k++)
	{
		uint offseti = placei - get<0>(anscentres[k]);
		uint offsetj = placej - get<1>(anscentres[k]);
		o=1;
			
			for (i=0; i<ansborders[k].size(); i++)
			{	
				uint new_x, new_y;
				tie(new_x, new_y) = ansborders[k][i];
				if ( (offseti + new_x > 0) && (offseti + new_x < out.n_rows) && (offsetj + new_y < out.n_cols) && (offsetj + new_y > 0)) {
					if (get<0>(out( offseti + new_x, offsetj + new_y) ) != 0) o=0;
				}
			}
			
		if (o==1) {Is_placed.push_back(1);} else {Is_placed.push_back(0);}	
	}
	
	uint ar=0;
	
	for (k=0; k<3; k++)
	{
		if ( (Is_placed[k]!=0) && (MAR[k]>ar) ) {ar=MAR[k]; result_idx=k;}
	}

// making final output image	
	uint offseti = placei - get<0>(anscentres[result_idx]);
	uint offsetj = placej - get<1>(anscentres[result_idx]);
	
	for (i=0; i<out.n_rows; ++i) {					
		for (j=0; j<out.n_cols; ++j)
		{
			if (get<0>(out(i,j))!=0)
			{
				get<1>(out(i,j)) = ((get<0>(out(i, j))*31)%255)+1;
				get<2>(out(i,j)) = ((get<0>(out(i, j))*7)%255)+1;
			}
		}
	}

	for (i=0; i<answers[result_idx].n_rows; ++i) {					
			for (j=0; j<answers[result_idx].n_cols; ++j)
			{
				
				if (get<0>( answers[result_idx](i,j) ) != 0)
				{
					tie(r, g, b) = answers[result_idx](i,j);
					r=255;
					g=0;
					b=0;
					out(i + offseti, j+offsetj) = make_tuple(r, g, b);
				}
				
			}
	}
	
	result_idx += 1;
	
    // Base: return array of found objects and index of the correct gear
    // Bonus: return additional parameters of gears
    return make_tuple(result_idx, object_array, out);
}

int main(int argc, char **argv)
{	
    if (argc != 4)
    {
        cout << "Usage: " << endl << argv[0]
             << " <in_image.bmp> <out_image.bmp> <out_result.txt>" << endl;
        return 0;
    }

    try {
        Image src_image = load_image(argv[1]);
        
        vector<Image> answers;
        string s1, s2, s3;        	
                
        s1 = argv[1];
        s2 = argv[1];
        s3 = argv[1];

        s1.insert(s1.find(".bmp"), "_1");
        s2.insert(s2.find(".bmp"), "_2");
        s3.insert(s3.find(".bmp"), "_3");

        const char* im1 = s1.c_str();
        const char* im2 = s2.c_str();
        const char* im3 = s3.c_str();

        Image img1 = load_image(im1);
		Image img2 = load_image(im2);
		Image img3 = load_image(im3);
		
		answers.push_back(img1);
		answers.push_back(img2);
		answers.push_back(img3);
        ofstream fout(argv[3]);

        vector<shared_ptr<IObject>> object_array;
        Image dst_image;
        int result_idx;
        tie(result_idx, object_array, dst_image) = repair_mechanism(src_image, answers);
        save_image(dst_image, argv[2]);

        fout << result_idx << endl;
        fout << object_array.size() << endl;
        for (const auto &obj : object_array)
            obj->Write(fout);

    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        return 1;
    }
}
