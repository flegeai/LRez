/*****************************************************************************
 *   LRez - Playing with linked reads
 *   
 *   BamComparator : Counts the common barcodes of regions from bam files
 *   Copyright (C) 2020  INRAE
 *   Authors: F. Legeai
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <stdexcept>
#include "hts.h"
#include "sam.h"
#include <map>
#include <list>
#include <vector>
#include <getopt.h>
#include <string.h>
//#include "boost/program_options.hpp"
#define BARCODE_SIZE 16
#define BXTAG "BX"

using namespace std;
//namespace po = boost::program_options;

// functions declarations
bool *fromflagtobits (int);
std::map<string,bool> getBarcodesfromRegion(samFile *samFile,  hts_idx_t  *index, bam_hdr_t *header, string region);
int common_barcodes(std::map<string,bool> list1, std::map<string,bool> list2 );

int main (int argc, char* argv[])
{

	try{
		string bam;
		string reg_files = "";
		string input = "";
		u_int size=1000;

		const char* const short_opts = "b:l:i:s:";

    static struct option long_options[] = {
        {"bam",      required_argument,       nullptr,  'b' },
        {"list", optional_argument,       nullptr,  'l' },
				{"in", optional_argument, nullptr , 'i' },
				{"size", optional_argument, nullptr , 's' },
        {nullptr,no_argument,nullptr,  0   }
    };

		while (true) {
			const auto opt = getopt_long(argc, argv,short_opts, long_options, nullptr);
			if (opt == -1) {break;}
		   switch (opt) {
             case 'b' :
						  	bam= std::string(optarg);
                 break;
						case 'l' :
							 reg_files= std::string(optarg);
								break;
						case 'i' :
								input=std::string(optarg);
								break;
						case 's' :
								size=	std::stoi(optarg);
								cerr << "Size set to: " << size << std::endl;
								break;

            default:
							std::cout <<
            	"--bam <bam>: bam file\n"
            	"--list <list>: list of regions (optional)\n"
							"--in : <string> name of a sequence \n"
							"--size <int> size of boundaries\n ";
                 exit(EXIT_FAILURE);
        }
    }

		if(!bam.empty()) {
			//open BAM for reading
			samFile *in = sam_open(bam.c_str(), "r");
			if(in == NULL) {
				throw runtime_error("Unable to open BAM/SAM file.");
			}
			//Load the index
			hts_idx_t *idx = sam_index_load(in, bam.c_str());
			if(idx == NULL) {
				throw runtime_error("Unable to open BAM/SAM index."
				" Make sure alignments are indexed");
			}

			//Get the header
			bam_hdr_t *header = sam_hdr_read(in);

			std::map<std::string,std::map<string,bool>> barcodes_map;

			// We extract barcodes from all the regions
			vector<string> regions;

			if (reg_files.compare("") != 0) {
				ifstream list_regions (reg_files);
				if (list_regions.is_open()){
					string region;
					while ( getline (list_regions, region) ){
						cerr << "Analyzing region "<< region << "\n";
						// mettre une map ici
						std::map<std::string, bool> l=getBarcodesfromRegion(in,idx,header,region);
						barcodes_map[region]=l;
						regions.push_back(region);
					}
				}

				// We generate the matrix
				cerr << "Producing the matrix\n";
				for(vector<string>::iterator r1= regions.begin(); r1 != regions.end();++r1){
					for(vector<string>::iterator r2= r1; r2 != regions.end();++r2){
						cout << *r1 << " " << *r2 << " "<< common_barcodes(barcodes_map[*r1],barcodes_map[*r2]) << "\n";
					}
				}
			}

			else {

				if (input.compare("") != 0) {

					std::map<std::string,std::map<string,bool>> barcodes_map;
					vector<string> regions;

					std::map<string,u_int> seq_length;
					// get the size of the segments
					for (int i = 0; i < header->n_targets; i++) {
						string region=header->target_name[i];
						if (header->target_len[i] < size) {
							string end = to_string(header->target_len[i]);
							region += ":0-" + end;
							//cout << region << "\n";
							std::map<std::string, bool> l=getBarcodesfromRegion(in,idx,header,region);
							barcodes_map[region]=l;
							regions.push_back(region);
						}
						else {
							// left
							string region_l = region + ":0-" + to_string(size);
							//cout << "left " << region_l << "\n";
							std::map<std::string, bool> l=getBarcodesfromRegion(in,idx,header,region_l);
							barcodes_map[region_l]=l;
							regions.push_back(region_l);
							// right
							u_int b_end = header->target_len[i]-size;
							string region_r = region + ":" + to_string(b_end) + "-" + to_string(header->target_len[i]);
							//cout << "right " << region_r << "\n";
						 	l=getBarcodesfromRegion(in,idx,header,region_r);
							barcodes_map[region_r]=l;
							regions.push_back(region_r);
						}
						seq_length[header->target_name[i]]=header->target_len[i];
					}

					string input_region = input;
					if (seq_length[input] < size) {
						input_region += ":0-" + to_string(seq_length[input]);
						for(vector<string>::iterator r= regions.begin(); r != regions.end();++r){
								int cb = common_barcodes(barcodes_map[input_region],barcodes_map[*r]);
								if (cb  > 0) {
									cout << input_region << " " << *r << " "<<  cb << "\n";
								}
						}
					}
					else {
							// left
							string input_region_l = input_region + ":0-" + to_string(size);
							for(vector<string>::iterator r= regions.begin(); r != regions.end();++r){
									int cb = common_barcodes(barcodes_map[input_region_l],barcodes_map[*r]);
									if (cb > 0 ) {
										cout << input_region_l << " " << *r << " "<< cb << "\n";
									}
							}
							// right
							u_int b_end = seq_length[input]-size;
							string input_region_r = input_region + ":" + to_string(b_end) + "-" + to_string(seq_length[input]);
							for(vector<string>::iterator r= regions.begin(); r != regions.end();++r){
									int cb = common_barcodes(barcodes_map[input_region_r],barcodes_map[*r]);
									if (cb > 0) {
										cout << input_region_r << " " << *r << " "<< cb << "\n";
									}
							}
					}
				}
			}
			hts_idx_destroy(idx);
			bam_hdr_destroy(header);
			sam_close(in);
				// Generation de	 la matrice

				return 0;
			}


	} catch (const runtime_error& e) {
		cerr << e.what();
	}
}


std::map<string,bool> getBarcodesfromRegion(samFile *samFile,  hts_idx_t  *index, bam_hdr_t *header, string region) {
	//Initialize iterator
	hts_itr_t *iter = NULL;
	//Move the iterator to the region we are interested in
	iter  = sam_itr_querys(index, header, region.c_str());
	if(header == NULL || iter == NULL) {
		sam_close(samFile);
		throw runtime_error("Unable to iterate to region within BAM.");
	}

	std::map<std::string,bool> barcode_map;
	//Initiate the alignment record
	bam1_t *aln = bam_init1();
	while(sam_itr_next(samFile, iter, aln) >= 0) {
		//printf("qname: %s\n",bam_get_qname(aln));
		//cout << "\tPos: " << aln->core.pos;
		//cout << "\tFlag: " << aln->core.flag;
		uint8_t* bxtag = bam_aux_get(aln, BXTAG);

		if (bxtag != 0) { // The  barcode is in the bam auxiliary tags (LongRanger)
			string barcode;
			int i = 1;

			while (*(bxtag+i) != '\0') {
				barcode+=*(bxtag+i);
				i++;
			}
			//cout<< barcode <<endl;
			barcode_map[barcode]=true;
		}
		else { // the barcode is in the sequence

			bool * bits = fromflagtobits(aln->core.flag);


			/* Barcode Extraction */
			/* We need to find the first seq of the fragments (bit 7 == true) */
			if (*(bits + 6) == true) {
				string barcode,seq;
				uint8_t *seqi = bam_get_seq(aln);
				if (*(bits +4) == false) { // Not reverse complemented
					for (int i = 0; i < BARCODE_SIZE; i++) {
						barcode += seq_nt16_str[bam_seqi(seqi, i)];
					}
				}
				else {
					for (int i = aln->core.l_qseq-1; i >= aln->core.l_qseq-BARCODE_SIZE; i--) {
						switch (seq_nt16_str[bam_seqi(seqi, i)]) {
							case 'T' :
							barcode+='A';
							break;
							case 'A' :
							barcode+='T';
							break;
							case 'G' :
							barcode+='C';
							break;
							case 'C' :
							barcode+='G';
							break;
							default :
							barcode+='N';
						}
					}
				}

				for (int i = 0; i < aln->core.l_qseq; i++) {
					seq += seq_nt16_str[bam_seqi(seqi, i)];
				}
				//cout << "\tSeq: " << seq << endl;
				//cout << "\tBarcode: " << barcode << endl;
				barcode_map[barcode]=true;
			}

		}
  	}
	// //	std::map<string,bool> barcodes;
	//
	// 	map<string, bool>::iterator itr;
	// 	for (itr = barcode_map.begin(); itr != barcode_map.end(); ++itr) {
	// 		//	cout <<  region  << " " << itr->first << "\n";
	// 			barcodes.push_back(itr->first);
	// 	}
	//
		hts_itr_destroy(iter);
		bam_destroy1(aln);
		return barcode_map;
}


bool* fromflagtobits(int n)
{
    static bool bits[12] ;
    int remainder;
		int step=0;
   while (n!=0)
    {
      	remainder = n%2;
      // 	cout << "Step " << step << ": " << n << "/2, Remainder = " << remainder << ", Quotient = " << n/2 << endl;
        n /= 2;
				bits[step]=remainder;
				step++;
    }
    return bits;
}

int common_barcodes(std::map<string,bool> list1, std::map<string,bool> list2 ) {
	int common =0;
	 for ( const auto &myPair : list1 ) {
		 common+=list2.count(myPair.first);
    //    std::cout << myPair.first << "\n";
    }

	// for (std::map<string,bool>::iterator it1 = list1.begin(); it1 != list1.end(); ++it1) {
	//
	// }
	return common;
}
