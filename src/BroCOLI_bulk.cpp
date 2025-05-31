#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <list>
#include <numeric>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <getopt.h>
#include <sys/stat.h>   
#include <sys/types.h>  
#include <dirent.h>
#include <thread>
#include <mutex>
#include <iomanip>
#include <condition_variable>
#include <queue>
#include <functional>
#include <atomic>
#include "thread_pool_executor.hpp"

std::mutex updatedGtfMutex;
std::mutex isoformCountMutex;
std::mutex geneCountMutex;
std::mutex traceMutex;
std::mutex bigMutex;


struct Split_Result{
    std::string read_name;
    std::vector<std::string> tokens;
    int read_length;
};


Split_Result get_string_split(const std::string& s, char delimiter) {
    Split_Result read_result = {};
    std::istringstream iss(s);
    std::string token;

    int number = 0;
    while (std::getline(iss, token, delimiter)) {
        number++;

        if (number == 1){
            read_result.read_name = token;
        } else if (1 < number && number <= 6) {

            read_result.tokens.push_back(token);
        } else if (number == 10) {
            read_result.read_length = token.size();
        }
    }
    return read_result;
}



int get_read_match_gene_length(std::vector<std::array<int,2>> CIGAR_vector){
    int Mlength = 0;
    for (auto each_array : CIGAR_vector) {
        Mlength = Mlength + (each_array[1] - each_array[0]);
    }
    return Mlength;
}


struct Read_intervals_and_Mlength{
    std::vector<std::array<int,2>> ReadIntervals;
    int ReadMatchLength;
};

Read_intervals_and_Mlength get_read_intervals(std::string CIGARvalue, std::string initpos){

    Read_intervals_and_Mlength CIGAR_Interval_Length;

    int MatchLength = 0;
    
    int position_begin = std::stoi(initpos);
    int position_end = 0;

    char currentChar = '0';
    int currentNumber = 0;
    bool inNumberMode = false;
    
    for (char ch : CIGARvalue) {

        if (std::isdigit(ch)) {
            inNumberMode = true;
            currentNumber = currentNumber * 10 + (ch - '0');
        } else {
            currentChar = ch;
            if (currentChar == 'M') {
                MatchLength = MatchLength + currentNumber;
                position_end = position_begin + currentNumber;
                std::array<int,2> each_small_QJ;
                each_small_QJ[0] = position_begin;
                each_small_QJ[1] = position_end;
                CIGAR_Interval_Length.ReadIntervals.push_back(each_small_QJ);

            } else if (currentChar == 'D') {

                position_end = position_begin + currentNumber;

            } else if (currentChar == 'I') {
                MatchLength = MatchLength + currentNumber;
                position_end = position_begin;

            } else if (currentChar == 'N') {
                position_end = position_begin + currentNumber;

            } else if (currentChar == 'S') {
                MatchLength = MatchLength + currentNumber;
                position_end = position_begin;

            } else if (currentChar == 'H') {
                MatchLength = MatchLength + currentNumber;
                position_end = position_begin;

            } else if (currentChar == 'P') {
                MatchLength = MatchLength + currentNumber;
                position_end = position_begin + currentNumber;

            } else {
                std::cout << currentChar << ' ' << "has not been considered yet!" << std::endl;
            }
            position_begin = position_end;
            inNumberMode = false;
            currentNumber = 0;
        }    
    }
    CIGAR_Interval_Length.ReadMatchLength = MatchLength;
    return CIGAR_Interval_Length;
}


struct Reads_Clusters {

    std::streampos lastPos;

    std::streampos newPos;

    std::map<std::string, std::vector<std::array<int,2>>> Mymap;

    std::map<std::string, int> Mylen;

    std::string SetRef_name;

    std::array<int,2> ClusterCoverage;

    int mapqless1 = 0;
    int mappingdiff = 0;
    int surveyNum = 0;    
}; 


Reads_Clusters get_each_cluster_reads(std::ifstream& samfile, std::streampos CurrentPos, std::streampos EndPos) {

    Reads_Clusters NewCluster = {};

    std::map<std::string, std::vector<std::array<int,2>>> read_informs;
    std::map<std::string, int> read_len;

    std::string line;
    std::string now_gene;
    std::string last_chr;

    samfile.seekg(CurrentPos, std::ios::beg);

    std::streampos earlyPos = CurrentPos;

    int early_begin_pos = 0;
    int early_end_pos = 0;

	while (getline(samfile, line))
	{

        earlyPos = CurrentPos; 
        CurrentPos = samfile.tellg();

        if (line[0] != '@'){
            NewCluster.surveyNum = NewCluster.surveyNum + 1;
            Split_Result This_Line = get_string_split(line, '\t');          

            int read_mapq = std::stoi(This_Line.tokens[3]) ;
            if (read_mapq > 1) {
                
                Read_intervals_and_Mlength CIGAR_interval;
                CIGAR_interval = get_read_intervals(This_Line.tokens[4], This_Line.tokens[2]);
                if (This_Line.read_length == CIGAR_interval.ReadMatchLength){

                    now_gene = This_Line.tokens[1];

                    int read_in_gene_begin_pos = std::stoi(This_Line.tokens[2]);
                    std::array<int,2> lastEle = CIGAR_interval.ReadIntervals.back();

                    int read_in_gene_end_pos = lastEle[1] - 1;

                    if (early_begin_pos == 0) {
                        early_begin_pos = read_in_gene_begin_pos;
                        early_end_pos = read_in_gene_end_pos;
                        last_chr = now_gene;
                    } else {

                        if (now_gene != last_chr){

                            NewCluster.lastPos = earlyPos;
                            NewCluster.newPos = CurrentPos;
                            NewCluster.Mymap = read_informs;
                            NewCluster.Mylen = read_len;
                            NewCluster.ClusterCoverage[0] = early_begin_pos;
                            NewCluster.ClusterCoverage[1] = early_end_pos;
                            NewCluster.SetRef_name = last_chr;
                            NewCluster.surveyNum = NewCluster.surveyNum - 1;
                            break;
                        } else{

                            if ((read_in_gene_begin_pos - early_end_pos > 0) || (early_begin_pos - read_in_gene_end_pos > 0)){

                                NewCluster.lastPos = earlyPos;
                                NewCluster.newPos = CurrentPos;

                                NewCluster.Mymap = read_informs;
                                NewCluster.Mylen = read_len;
                                NewCluster.ClusterCoverage[0] = early_begin_pos;
                                NewCluster.ClusterCoverage[1] = early_end_pos;
                                NewCluster.SetRef_name = last_chr;
                                NewCluster.surveyNum = NewCluster.surveyNum - 1;
                                break;
                            } else {

                                early_begin_pos = (read_in_gene_begin_pos>early_begin_pos)?early_begin_pos:read_in_gene_begin_pos;
                                early_end_pos = (read_in_gene_end_pos>early_end_pos)?read_in_gene_end_pos:early_end_pos;
                                last_chr = now_gene;
                            }
                        }
                    }
                    read_informs[This_Line.read_name] = CIGAR_interval.ReadIntervals;
                    read_len[This_Line.read_name] = This_Line.read_length;

                } else {
                    NewCluster.mappingdiff = NewCluster.mappingdiff + 1;
                }
                
            } else {
                NewCluster.mapqless1 = NewCluster.mapqless1 + 1;
            }
        }
	}

    if (CurrentPos >= EndPos) {

        NewCluster.lastPos = earlyPos;
        NewCluster.newPos = CurrentPos;

        NewCluster.Mymap = read_informs;
        NewCluster.Mylen = read_len;
        NewCluster.ClusterCoverage[0] = early_begin_pos;
        NewCluster.ClusterCoverage[1] = early_end_pos;
        NewCluster.SetRef_name = last_chr;
    }

    return NewCluster;
}



struct SpliceJs {

    std::unordered_map<std::string, std::vector<std::array<int,2>>> reads_SJs;

    std::unordered_map<std::string, std::vector<std::array<int,2>>> reads_SJs_left;

    std::unordered_map<std::string, std::vector<std::array<int,2>>> reads_SJs_right;

    std::unordered_map<std::string, std::array<int,2>> reads_begin_end;

    std::unordered_map<std::string, std::array<int,2>> reads_single_exon;

}; 


SpliceJs get_reads_allSJs(std::map<std::string, std::vector<std::array<int,2>>>& RCs, int SJD) {

    SpliceJs NewSJs = {};

    std::unordered_map<std::string, std::vector<std::array<int,2>>> AllRead_SJs;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> AllRead_SJs_left;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> AllRead_SJs_right;
    std::unordered_map<std::string, std::array<int,2>> AllRead_begin_end;
    std::unordered_map<std::string, std::array<int,2>> AllRead_SingleExon;

    for (const auto &pair : RCs) {

        std::string Rname = pair.first; 
        std::vector<std::array<int,2>> RMs = pair.second;

        std::vector<std::array<int,2>> ReadSJs;
        std::vector<std::array<int,2>> ReadSJs_left;
        std::vector<std::array<int,2>> ReadSJs_right;

        std::array<int,2> Read_SingleExon;
        std::array<int,2> Read_begin_end;

        for (auto it = RMs.begin(); it!= RMs.end()-1; it++) {

            if (it == RMs.begin()){
                Read_begin_end[0] = (*(it))[0];
            }
            if (it== RMs.end()-2){
                Read_begin_end[1] = (*(it+1))[1];
            }
            int later = (*(it+1))[0];
            int former = (*(it))[1];
            if (later - former > SJD) {

                std::array<int,2> eachSJ;
                eachSJ[0] = former;
                eachSJ[1] = later - 1;
                ReadSJs.push_back(eachSJ);

                std::array<int,2> eachSJ_left;
                eachSJ_left[0] = (*(it))[0];
                eachSJ_left[1] = former - 1;
                ReadSJs_left.push_back(eachSJ_left);
  
                std::array<int,2> eachSJ_right;
                eachSJ_right[0] = later;
                eachSJ_right[1] = (*(it+1))[1];
                ReadSJs_right.push_back(eachSJ_right);
            }
        }

        if (ReadSJs.size() == 0) {
            Read_SingleExon[0] = (*RMs.begin())[0];
            Read_SingleExon[1] = (*(RMs.end()-1))[1] - 1;
            AllRead_SingleExon[Rname] = Read_SingleExon;
        } else {

            AllRead_SJs[Rname] = ReadSJs;
            AllRead_SJs_left[Rname] = ReadSJs_left;
            AllRead_SJs_right[Rname] = ReadSJs_right;
            AllRead_begin_end[Rname] = Read_begin_end;
        }

    }

    NewSJs.reads_SJs = AllRead_SJs;
    NewSJs.reads_SJs_left = AllRead_SJs_left;
    NewSJs.reads_SJs_right = AllRead_SJs_right;
    NewSJs.reads_begin_end = AllRead_begin_end;
    NewSJs.reads_single_exon = AllRead_SingleExon;
    return NewSJs;
}



struct GTF
{
    std::map<std::string, std::map<std::string, std::vector<std::array<int,2>>>> GTF_transcript;
    std::map<std::string, std::map<std::string, std::array<int,2>>> GTF_gene; 
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> GTF_transcript_strand;
    std::map<std::string, std::map<std::string, std::string>> GTF_gene_strand;
    std::map<std::string, std::map<std::string, std::vector<std::string>>> GTF_gene2transcript;
};

GTF get_gtf_annotation(std::string& GTFFile_name){
    GTF GTFAll_Info;
    if (!GTFFile_name.empty()){

        std::cout << "***** Now open the gtf file: " << GTFFile_name << "! *****" << std::endl;        

        std::map<std::string, std::vector<std::array<int,2>>> ChrEach;
        std::unordered_map<std::string, std::string> ChrTranscriptStrand;
        std::map<std::string, std::array<int,2>> GeneEach;      
        std::vector<std::array<int,2>> GTFAnno_SJs;  
        std::array<int,2> annoExon;
        
        std::ifstream GTFFile; 

        GTFFile.open(GTFFile_name);
        if (!GTFFile.is_open())	{
            std::cout << "We can't open GTF file !" << GTFFile_name << std::endl;
            exit(EXIT_FAILURE);
        }        
        std::string line;
        std::string now_chr_name = "Han";
        std::string early_chr_name = "Han";
        std::string early_gene_transcript_name = "WHH";
        std::string now_gene_transcript_name = "WHH";

        std::string early_Exonstrand = "Wei";
        std::string now_Exonstrand = "Wei";

        while (std::getline(GTFFile, line)){
          
            if (line[0] != '#'){

                std::istringstream lineiss(line);

                std::vector<std::string> tokens;

                std::string token;

                while (std::getline(lineiss, token, '\t')) {
                    tokens.push_back(token);
                } 
                
                if (tokens[2] == "exon"){
                    early_chr_name = now_chr_name;
                    now_chr_name = tokens[0]; 
                    early_Exonstrand = now_Exonstrand;
                    now_Exonstrand = tokens[6]; 

                    std::string AllEndT = tokens.back();

                    std::istringstream EndTT(AllEndT);
                    std::string Endtoken;
                    std::vector<std::string> Endtokens;

                    while (std::getline(EndTT, Endtoken, '"')){
                        Endtokens.push_back(Endtoken);
                    }
                    early_gene_transcript_name = now_gene_transcript_name;
                    now_gene_transcript_name = Endtokens[1]+'|'+Endtokens[3];

                    if ((GTFAnno_SJs.size() != 0) && (early_gene_transcript_name != now_gene_transcript_name)){
                        if ((GTFAnno_SJs.size()>1) && (GTFAnno_SJs[0][0] > GTFAnno_SJs[1][0])){ 
                            std::reverse(GTFAnno_SJs.begin(),GTFAnno_SJs.end());
                        } 
                        ChrEach[early_gene_transcript_name] = GTFAnno_SJs;
                        ChrTranscriptStrand[early_gene_transcript_name] = early_Exonstrand;
                        GTFAnno_SJs.clear();                        
                    }

                    if ((now_chr_name != early_chr_name) && (ChrEach.size() != 0)){
                        GTFAll_Info.GTF_transcript[early_chr_name] = ChrEach;
                        GTFAll_Info.GTF_transcript_strand[early_chr_name] = ChrTranscriptStrand;
                        ChrEach.clear();
                        ChrTranscriptStrand.clear();                     
                    }
                    annoExon[0] = std::stoi(tokens[3]);
                    annoExon[1] = std::stoi(tokens[4]);
                    GTFAnno_SJs.push_back(annoExon);
                }
            }
        }

        if (GTFFile.eof()) {

            if ((GTFAnno_SJs.size()>1) && (GTFAnno_SJs[0][0] > GTFAnno_SJs[1][0])) {
                std::reverse(GTFAnno_SJs.begin(), GTFAnno_SJs.end());
            }            
            ChrEach[now_gene_transcript_name] = GTFAnno_SJs;
            GTFAll_Info.GTF_transcript[now_chr_name] = ChrEach;
            ChrTranscriptStrand[now_gene_transcript_name] = now_Exonstrand;
            GTFAll_Info.GTF_transcript_strand[now_chr_name] = ChrTranscriptStrand;
            std::cout << "***** The GTF file has been read to the end ! *****" << std::endl;
        } 
        else if (GTFFile.fail()) {
            std::cout << "File FALSE !" << std::endl;
        }
        else {std::cout << "unkown reason" << std::endl;}

        GTFFile.close();

        std::vector<int> All_gene_begin;
        std::vector<int> All_gene_end;
        std::vector<std::string> tx_name;
        size_t pos;
        std::string now_gene_name;
        std::string ago_gene_name;
        std::string now_chr;
        std::string ago_chr;
        std::string Strand;
     
        for (const auto& eachChr:GTFAll_Info.GTF_transcript){

            now_chr = eachChr.first;
            GTFAll_Info.GTF_gene[now_chr] = {};
            GTFAll_Info.GTF_gene_strand[now_chr] = {};
            GTFAll_Info.GTF_gene2transcript[now_chr] = {};

            if (All_gene_begin.size() != 0 && now_chr != ago_chr) {
                auto min_it = std::min_element(All_gene_begin.begin(), All_gene_begin.end());
                auto max_it = std::max_element(All_gene_end.begin(), All_gene_end.end());
                GTFAll_Info.GTF_gene[ago_chr][ago_gene_name][0] = *min_it;
                GTFAll_Info.GTF_gene[ago_chr][ago_gene_name][1] = *max_it;
                GTFAll_Info.GTF_gene2transcript[ago_chr][ago_gene_name] = tx_name;
                All_gene_begin.clear();
                All_gene_end.clear();
                tx_name.clear();
                ago_gene_name = "";
            }

            for (const auto& eachGene:eachChr.second) {
                size_t pos = eachGene.first.find('|');
                now_gene_name = eachGene.first.substr(0, pos);
                Strand = GTFAll_Info.GTF_transcript_strand[now_chr][eachGene.first];

                if (ago_gene_name == now_gene_name) {
                    All_gene_begin.push_back(eachGene.second[0][0]);
                    All_gene_end.push_back(eachGene.second[eachGene.second.size()-1][1]);
                    tx_name.push_back(eachGene.first);
                } else {
                    if (ago_gene_name.size() != 0) {
                        auto min_it = std::min_element(All_gene_begin.begin(), All_gene_begin.end());
                        auto max_it = std::max_element(All_gene_end.begin(), All_gene_end.end());
                        GTFAll_Info.GTF_gene[eachChr.first][ago_gene_name][0] = *min_it;
                        GTFAll_Info.GTF_gene[eachChr.first][ago_gene_name][1] = *max_it;
                        GTFAll_Info.GTF_gene2transcript[eachChr.first][ago_gene_name] = tx_name;
                    }
                    All_gene_begin.clear();
                    All_gene_end.clear();
                    tx_name.clear();
                    All_gene_begin.push_back(eachGene.second[0][0]);
                    All_gene_end.push_back(eachGene.second[eachGene.second.size()-1][1]);
                    tx_name.push_back(eachGene.first);
                    GTFAll_Info.GTF_gene_strand[now_chr][now_gene_name] = Strand;
                }
                ago_gene_name = now_gene_name;
            }
            ago_chr = now_chr;
        }
    } else {
        std::cout << "***** No GTF files are put into the program! *****" << std::endl;         
    }
    return GTFAll_Info;
}



struct GTFsj{

    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::array<int,2>>>> mSJs;

    std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> mSJsBE;

    std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> SE;
};


GTFsj get_SJs_SE(std::map<std::string, std::map<std::string, std::vector<std::array<int,2>>>>& Known_Exon){

    GTFsj Known_SJ_SE;

    if (Known_Exon.size() != 0){

        std::map<std::string, std::vector<std::array<int,2>>> Every_transcript_All;
        std::vector<std::array<int,2>> Every_transcript;
        std::string transaction_name;
        std::string chr_name;
        
        std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::array<int,2>>>> chr_Every_SJs_all;
        std::unordered_map<std::string, std::vector<std::array<int,2>>> Every_SJs_All;
        std::vector<std::array<int,2>> Every_SJs;
        std::array<int,2> SSJs;

        std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> chr_Every_SE_all;
        std::unordered_map<std::string, std::array<int,2>> Every_SE_All;
        std::array<int,2> Every_SE;

        std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> chr_Every_SJ_begin_end_all;
        std::unordered_map<std::string, std::array<int,2>> Every_SJ_begin_end_all;
        std::array<int,2> Every_SJ_begin_end;

        for (const auto& pair : Known_Exon) {

            chr_name = pair.first;
            Every_transcript_All = pair.second;

            Every_SJs_All.clear();
            Every_SE_All.clear();
            Every_SJ_begin_end_all.clear();

            for (auto it = Every_transcript_All.begin(); it != Every_transcript_All.end(); ++it) {

                transaction_name = it->first;
                Every_transcript = it->second;

                if (Every_transcript.size() == 1){

                    Every_SE[0] = (Every_transcript.front())[0];
                    Every_SE[1] = (Every_transcript.front())[1];
                    Every_SE_All[transaction_name] = Every_SE;
                } else {

                    Every_SJs.clear();
                    for (auto exonit = Every_transcript.begin(); exonit != Every_transcript.end() - 1; ++exonit) {
  
                        if (exonit == Every_transcript.begin()) {
                            Every_SJ_begin_end[0] = (*exonit)[0];
                        }
                        if (exonit == Every_transcript.end() - 2) {
                            Every_SJ_begin_end[1] = ((*(exonit+1))[1]);
                        }
                        SSJs[0] = ((*exonit)[1])+1;
                        SSJs[1] = ((*(exonit+1))[0])-1;
                        Every_SJs.push_back(SSJs);
                    }
                    Every_SJs_All[transaction_name] = Every_SJs;
                    Every_SJ_begin_end_all[transaction_name] = Every_SJ_begin_end;
                }
            }

            if (Every_SJs_All.size() != 0) {
                chr_Every_SJs_all[chr_name] = Every_SJs_All;
            }
            if (Every_SE_All.size() != 0) {
                chr_Every_SE_all[chr_name] = Every_SE_All;
            }
            if (Every_SJ_begin_end_all.size() != 0) {
                chr_Every_SJ_begin_end_all[chr_name] = Every_SJ_begin_end_all;
            }
        }

        Known_SJ_SE.mSJs = chr_Every_SJs_all;
        Known_SJ_SE.SE = chr_Every_SE_all;
        Known_SJ_SE.mSJsBE = chr_Every_SJ_begin_end_all;
    }
    return Known_SJ_SE;
}



std::unordered_map<std::string, std::string> Read_fasta_file(std::string& fastafile_name){

    std::unordered_map<std::string, std::string> chrgeneString;

    std::vector<std::string> VecStrings;

    std::cout << "***** Now open Fasta file: " << fastafile_name << "! *****" << std::endl;

    std::ifstream FastaFile; 

	FastaFile.open(fastafile_name);

	if (!FastaFile.is_open()) {
		std::cout << "We can't open Fasta file !" << fastafile_name << std::endl;
		exit(EXIT_FAILURE);
	}

    std::string ChrGeneName;
    std::string concatenatedString;
    std::stringstream accumulatedStringStream;
    std::string line;
    std::string token;
    
    while (std::getline(FastaFile, line))
    {      
        if (line[0] == '>'){

            if (!ChrGeneName.empty()){
                concatenatedString = accumulatedStringStream.str();
                chrgeneString[ChrGeneName] = concatenatedString;
                accumulatedStringStream.str("");
            }
            std::istringstream iss(line);
            int count = 0;
            while (std::getline(iss, token, ' ')){
                count++;

                if (count == 1){
                    ChrGeneName = token.substr(1);
                }
            }

        } else {

            accumulatedStringStream << line;
        }
    }

    if (!ChrGeneName.empty()) {

        concatenatedString = accumulatedStringStream.str();
        chrgeneString[ChrGeneName] = concatenatedString;
    }

    if (FastaFile.eof()) {

        std::cout << "***** The Fasta file has been read to the end ! *****" << std::endl;
    } 
    else if (FastaFile.fail()) {
        std::cout << "File FALSE !" << std::endl;
    }
    else {std::cout << "unkown reason" << std::endl;}

    FastaFile.close();
    return chrgeneString;
}



bool ifSjNoError(const std::vector<std::array<int,2>>& leftSj, const std::vector<std::array<int,2>>& rightSj, int& Sj_Number) {
    bool flag;
    std::array<int,2> thisSj_left = leftSj[Sj_Number];
    std::array<int,2> thisSj_right = rightSj[Sj_Number];
    if ((thisSj_left[1]-thisSj_left[0]>10) && (thisSj_right[1]-thisSj_right[0]>10)) {
        flag = true;
    } else {
        flag = false;
    }
    return flag;
}


int ifSjSignal(std::string& Fasta_chr, const std::array<int,2>& Sj_Array) {
    int flag = 0;
    std::string signal1 = Fasta_chr.substr(Sj_Array[0]-1, 2);
    std::string signal2 = Fasta_chr.substr(Sj_Array[1]-2, 2);

    std::transform(signal1.begin(), signal1.end(), signal1.begin(), ::toupper);
    std::transform(signal2.begin(), signal2.end(), signal2.begin(), ::toupper);

    if (signal1 == "GT" && signal2 == "AG"){
        flag = 1;

    } else if (signal1 == "CT" && signal2 == "AC"){
        flag = 2;

    } else if (signal1 == "GC" && signal2 == "AG"){
        flag = 1;

    } else if (signal1 == "CT" && signal2 == "GC"){
        flag = 2;

    } else if (signal1 == "AT" && signal2 == "AC"){
        flag = 1;

    } else if (signal1 == "GT" && signal2 == "AT"){
        flag = 2;

    } else {
        flag = 0;
    }                
    return flag;
}



static struct option long_options[] = {
    {"sam", required_argument, 0, 's'},
    {"fasta", required_argument, 0, 'f'},
    {"gtf", required_argument, 0, 'g'},
    {"output", required_argument, 0, 'o'},
    {"mode", required_argument, 0, 'm'},
    {"SJDistance", required_argument, 0, 'j'},  
    {"support", required_argument, 0, 'n'},   
    {"single_exon_boundary", required_argument, 0, 'e'},
    {"graph_distance", required_argument, 0, 'd'},  
    {"thread", required_argument, 0, 't'}, 
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}  
};

bool directoryExists(const std::string& path) {
    struct stat info;

    if (stat(path.c_str(), &info) != 0) {
        return false;
    } else if (info.st_mode & S_IFDIR) {
        return true;
    } else {
        return false;
    }
}

bool createDirectory(const std::string& path) {
    if (mkdir(path.c_str(), 0755) == 0) {
        return true;
    } else {
        std::cerr << "Error creating directory: " << path << strerror(errno) << std::endl;
        return false;
    }
}

std::string joinPath(const std::string& directory, const std::string& filename) {
    if (directory.back() == '/') {
        return directory + filename;  
    } else {
        return directory + "/" + filename;  
    }
}


void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << std::endl;
    std::cout << "Required parameter: [-s <sam_file>] [-f <fasta_file>] [-o <output_file>] " << std::endl;
    std::cout << "Optional parameters: [-g <gtf_file>] [-j <SJDistance>] [-n <SJ_support_read_number>] [-d <Graph_distance>] [-t <Thread>] [-m <mode>]" << std::endl;
    std::cout << "Help: [-h <help>]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -s, --sam                   SAM file path. We recommend using absolute paths. If you have a single file, you can directly provide its absolute path. If you have multiple files, you can specify the path to a folder that contains all the sorted SAM files you want to process. (required)" << std::endl;
    std::cout << "  -f, --fasta                 FASTA file path. FASTA file requires the chromosome names to match the GTF file. (required)" << std::endl;
    std::cout << "  -g, --gtf                   input annotation file in GTF format. (optional, Recommendation provided)" << std::endl;
    std::cout << "  -o, --output                output folder path. (required)" << std::endl;
    std::cout << "  -m, --mode                  Sequencing method (default: 0(cDNA), !0(direct RNA))" << std::endl;
    std::cout << "  -j, --SJDistance            the minimum distance determined as intron. (optional, default:18)" << std::endl;
    std::cout << "  -n, --support               min perfect read count for all splice junctions of novel isoform. (optional, default:2)" << std::endl;
    std::cout << "  -e, --single_exon_boundary  belongs to the isoform scope of a single exon. (optional, default:60)" << std::endl;
    std::cout << "  -d, --graph_distance        the distance threshold for constructing the isoform candidate distance graph. (optional, default:60)" << std::endl;
    std::cout << "  -t, --thread                thread number. (optional, default:8)" << std::endl;
    std::cout << "  -h, --help                  show this help information." << std::endl;
}

std::vector<std::string> check_catalog_exist(const std::string& output_path,
                                            std::vector<std::string>& samFileVec) {
    std::vector<std::string> outputFileVector;

    if (!directoryExists(output_path)) {
        std::cout << "Directory does not exist, creating it..." << std::endl;

        if (createDirectory(output_path)) {
            std::cout << "Directory created successfully: " << output_path << std::endl;
        } else {
            std::cerr << "Failed to create directory: " << output_path << std::endl;
            exit(EXIT_FAILURE);
        }
    } else {
        std::cout << "Directory already exists: " << output_path << std::endl;
    } 


    std::string updatedGtfPath = joinPath(output_path, "updated_annotitions.gtf");
    outputFileVector.push_back(updatedGtfPath);
    std::ofstream gtf_file(updatedGtfPath, std::ios::trunc);
    gtf_file.close();

    std::string quantTranscriptPath = joinPath(output_path, "counts_transcript.txt");
    outputFileVector.push_back(quantTranscriptPath);
    std::ofstream Q_output_Transcript(quantTranscriptPath, std::ios::trunc);
    
    std::string quantGenePath = joinPath(output_path, "counts_gene.txt");
    outputFileVector.push_back(quantGenePath);
    std::ofstream Q_output_Gene(quantGenePath, std::ios::trunc);
    
    std::string TracePath = joinPath(output_path, "compatible_isoform.tsv");
    outputFileVector.push_back(TracePath);
    std::ofstream TraceIsoform(TracePath, std::ios::trunc);

    if (samFileVec.size() == 1){
        if (Q_output_Transcript.is_open()){
            Q_output_Transcript << "transcript_id" << '\t' << "gene_id" << '\t' << samFileVec[0] << '\n';
        }
        if (Q_output_Gene.is_open()){
            Q_output_Gene << "gene_id" << '\t' << samFileVec[0] << '\n';
        }
        if (TraceIsoform.is_open()){
            TraceIsoform << "read_id" << '\t' << "category" << '\t' << "isoform_id" << '\t' << "gene_id" << '\t' << "file" << '\n';
        }
    } else {
        if (Q_output_Transcript.is_open()){
            Q_output_Transcript << "transcript_id" << '\t' << "gene_id";
        }
        if (Q_output_Gene.is_open()){
            Q_output_Gene << "gene_id";
        }        
        if (TraceIsoform.is_open()){
            TraceIsoform << "read_id" << '\t' << "category" << '\t' << "isoform_id" << '\t' << "gene_id" << '\t' << "file" << '\n';
        }
        for (const auto& eachFile:samFileVec) {
            Q_output_Transcript << '\t' << eachFile;
            Q_output_Gene << '\t' << eachFile;
        }
        Q_output_Transcript << '\n';
        Q_output_Gene << '\n';
    }
    Q_output_Transcript.close();
    Q_output_Gene.close(); 
    TraceIsoform.close();   
    return outputFileVector;
}



std::vector<std::string> traverse_sam_file(const std::string& sam_file_path, const std::string& output_path){

    std::vector<std::string> sam_file_vector;

    struct stat sam_stat;
 
    stat(sam_file_path.c_str(), &sam_stat);

    if (S_ISREG(sam_stat.st_mode)) {

        std::cout << "* Only one sam file is entered! * << " << sam_file_path << std::endl;

        sam_file_vector.push_back(sam_file_path);

    } else if (S_ISDIR(sam_stat.st_mode)) {

        std::cout << "* A folder was entered! * << " << sam_file_path << std::endl;
        DIR* dir = opendir(sam_file_path.c_str());
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            if (entry->d_name[0] != '.') { 
                std::string fileName = entry->d_name;
                std::cout << "Find file: " << fileName << std::endl;

                if (fileName.rfind(".sam") == fileName.length() - 4) {
                    sam_file_vector.push_back(fileName);
                }
            }
        }
        closedir(dir);
        std::cout << "^-^ There are " << sam_file_vector.size() << " sam files in total. ^-^" << std::endl;
    } else {
        std::cerr << "* Not a valid file or folder! *" << strerror(errno) << std::endl;
    }
    std::string File_explain_path = joinPath(output_path, "file_explain.txt");
    std::ofstream Explain_file(File_explain_path, std::ios::trunc);    
    Explain_file << "File" << '\t' << "File_Path" << '\n';
    for (int i = 0; i < sam_file_vector.size(); i++) {
        std::cout << "SAM File " << i << " : " << sam_file_vector[i] << std::endl;
        Explain_file << i << '\t' << sam_file_vector[i] << '\n';
    }
    Explain_file.close();
    return sam_file_vector;
}


std::streampos findNextLineStart(std::ifstream& file, std::streampos pos) {
    file.seekg(pos);
    std::string line;

    if (std::getline(file, line)) {
        return file.tellg(); 
    }
    return pos; 
}



void processChunk(const std::string& one_sam_file_path, const std::streampos& start, const std::streampos& end, const std::string& output_path, const int file_i, const int SJ_Distance, std::unordered_map<std::string, std::string>& FastaRef){

    int Group_index = 0;
    std::string chrchr;
    Reads_Clusters each_cluster_informs;
    SpliceJs each_read_SJs_informs;
    std::array<int,2> readBeginEnd;
    
    bool thisSjNoError = 0;
    int thisSjSignal = 0;
    int SjNumber = 0; 

    std::ostringstream oss;
    oss << "Read_" << std::setw(4) << std::setfill('0') << file_i << ".txt";
    std::string File_name = oss.str();
    std::string ReadInformPath = joinPath(output_path, File_name);
    std::ofstream ReadInform(ReadInformPath, std::ios::trunc);


    std::ifstream samfile(one_sam_file_path);
    samfile.seekg(start);

    std::streampos Current_Position = samfile.tellg();
    std::streampos Last_Position = Current_Position;
    
    while (Current_Position < end) {
        Group_index++;
        each_cluster_informs = get_each_cluster_reads(samfile, Last_Position, end);
        chrchr = each_cluster_informs.SetRef_name;
        Last_Position = each_cluster_informs.lastPos;
        Current_Position = each_cluster_informs.newPos;

        each_read_SJs_informs = get_reads_allSJs(each_cluster_informs.Mymap, SJ_Distance);  
                                                                                            
        for (const auto& eachRead:each_read_SJs_informs.reads_SJs) {
            readBeginEnd = each_read_SJs_informs.reads_begin_end[eachRead.first];
            ReadInform << eachRead.first << '\t' << chrchr << '\t' << Group_index << '\t'
                        << each_cluster_informs.ClusterCoverage[0] << '\t' 
                        << each_cluster_informs.ClusterCoverage[1] << '\t' << readBeginEnd[0] << '\t' 
                        << readBeginEnd[1] << '\t' << each_cluster_informs.Mylen[eachRead.first];
            
            SjNumber = 0;
            for (const auto& Sj:eachRead.second) {
                ReadInform << '\t' << Sj[0] << '\t' << Sj[1];

                thisSjNoError = ifSjNoError(each_read_SJs_informs.reads_SJs_left[eachRead.first], each_read_SJs_informs.reads_SJs_right[eachRead.first], SjNumber); 
                if (thisSjNoError == 1) {

                    thisSjSignal = ifSjSignal(FastaRef[chrchr], Sj);
                    if (thisSjSignal == 1) {
                        ReadInform << '\t' << 1;
                    } else if (thisSjSignal == 2) {
                        ReadInform << '\t' << 2;
                    } else {
                        ReadInform << '\t' << 0;
                    }
                } else {
                    ReadInform << '\t' << 0;
                }
                SjNumber++;                    
            }
            ReadInform << '\n';
        }

        for (const auto& eachRead:each_read_SJs_informs.reads_single_exon) {
            readBeginEnd = eachRead.second;
            ReadInform << eachRead.first << '\t' << chrchr << '\t' << Group_index << '\t'
                        << each_cluster_informs.ClusterCoverage[0] << '\t' 
                        << each_cluster_informs.ClusterCoverage[1] << '\t' << readBeginEnd[0] << '\t' 
                        << readBeginEnd[1] << '\t' << each_cluster_informs.Mylen[eachRead.first] << '\n';            
        } 
    } 
    samfile.close();
    ReadInform.close();
    std::cout << "^-^ Thread: " << std::this_thread::get_id() << " has completed processing! ^-^" << std::endl;
}



Split_Result get_line_split(const std::string& s, char delimiter){

    Split_Result read_result = {};

    std::istringstream iss(s);

    std::string token;

    int number = 0;
    while (std::getline(iss, token, delimiter)){
        number++;
        if (number == 1) {
            read_result.read_name = token;
        } else {
            read_result.tokens.push_back(token);
        }
    }
    return read_result;
}


Split_Result get_sj_split(const std::string& s, char delimiter){

    Split_Result sj_result = {};

    std::istringstream iss(s);

    std::string token;

    while (std::getline(iss, token, delimiter)){
        sj_result.tokens.push_back(token);
    }
    return sj_result;
}



struct FileSplit {

    std::map<std::string, std::vector<std::array<int,2>>> chr_coverage; 
    std::map<std::array<int,2>, std::array<std::streampos,2>> coverage2pos;

    std::vector<std::streampos> reads_pointer; 

    std::vector<int> group_reads_number; 

    std::string readtxt_path; 
    int FileNo;
};

template <typename T>
std::vector<std::size_t> sort_indexes_e(std::vector<T> &v)
{
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    return idx;
}



FileSplit Merge_Read_Small_Files(const std::string& SmallFilePath, const int& SAMFileNumber) {
    FileSplit Chunk_Bang;

    std::vector<std::string> Read_x_vec;

    DIR* dir = opendir(SmallFilePath.c_str());
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_name[0] == 'R') { 
            Read_x_vec.push_back(entry->d_name);
        } 
    }
    std::sort(Read_x_vec.begin(), Read_x_vec.end()); 
    closedir(dir);

    std::streampos Readpos;
    std::string File_total_name = SmallFilePath + "/All_Read.txt"; 
    std::ofstream AllReadInform(File_total_name, std::ios::trunc);  
    Chunk_Bang.readtxt_path = File_total_name;
    Readpos = AllReadInform.tellp();
    Chunk_Bang.reads_pointer.push_back(Readpos);

    std::string line;
    std::string chr_name;
    std::string last_chr;
    std::string earlyGroup = "+";
    std::string thisGroup = "+";
    int Group_index = 0;
    int Group_low;
    int Group_high;

    std::map<std::string, std::vector<std::string>> Group_temporary;
    std::array<int,2> this_coverage;

    for (int file_number = 0; file_number < Read_x_vec.size(); file_number++){
        std::string fileName = SmallFilePath + "/" + Read_x_vec[file_number];
        std::ifstream SmallSamFile(fileName);
        earlyGroup = "+";
        thisGroup = "+";

        while (getline(SmallSamFile, line)) {
            Split_Result This_Line = get_line_split(line, '\t');
            thisGroup = This_Line.tokens[1];
            this_coverage = {std::stoi(This_Line.tokens[2]), std::stoi(This_Line.tokens[3])};

            if (thisGroup == earlyGroup) {

                auto it = Group_temporary.find(This_Line.read_name);
                if (it != Group_temporary.end()){
                    if (std::stoi(This_Line.tokens[6]) > std::stoi(Group_temporary[This_Line.read_name][6])){
                        Group_temporary[This_Line.read_name] = This_Line.tokens;
                    }
                } else {
                    Group_temporary[This_Line.read_name] = This_Line.tokens;
                }
            } else { 

                if (Group_temporary.size() != 0) {

                    if (last_chr == This_Line.tokens[0]) {

                        if (Group_low <= this_coverage[1] && this_coverage[0] <= Group_high){

                            auto it = Group_temporary.find(This_Line.read_name);
                            if (it != Group_temporary.end()){
                                if (std::stoi(This_Line.tokens[6]) > std::stoi(Group_temporary[This_Line.read_name][6])){
                                    Group_temporary[This_Line.read_name] = This_Line.tokens;
                                }
                            } else {
                                Group_temporary[This_Line.read_name] = This_Line.tokens;
                            }
                            Group_high = this_coverage[1]>Group_high ? this_coverage[1] : Group_high;
                            Group_low = this_coverage[0]<Group_low ? this_coverage[0] : Group_low;
                        } else {

                            Group_index++;
                            for (const auto& Read:Group_temporary){
                                AllReadInform << SAMFileNumber << '\t' << Read.first;
                                for (int j = 0; j < Read.second.size(); j++) {
                                    if (j == 1) {
                                        AllReadInform << '\t' << Group_index;
                                    } else if (j == 2) {
                                        AllReadInform << '\t' << Group_low;
                                    } else if (j == 3) {
                                        AllReadInform << '\t' << Group_high;
                                    } else {
                                        AllReadInform << '\t' << Read.second[j];
                                    }
                                }
                                chr_name = Read.second[0];
                                AllReadInform << '\n';
                            }

                            Readpos = AllReadInform.tellp();
                            Chunk_Bang.coverage2pos.insert({{Group_low, Group_high}, {Chunk_Bang.reads_pointer.back(),Readpos}});
                            Chunk_Bang.reads_pointer.push_back(Readpos); 
                            Chunk_Bang.reads_pointer.push_back(Readpos); 
                            Chunk_Bang.group_reads_number.push_back(Group_temporary.size());

                            auto it = Chunk_Bang.chr_coverage.find(chr_name);
                            if (it != Chunk_Bang.chr_coverage.end()){
                                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                            } else {
                                Chunk_Bang.chr_coverage[chr_name] = {};
                                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                            }
                            
                            Group_temporary.clear();
                            Group_temporary[This_Line.read_name] = This_Line.tokens;
                            Group_high = this_coverage[1];
                            Group_low = this_coverage[0];                            
                        }
                    } else {

                        Group_index++;
                        for (const auto& Read:Group_temporary){
                            AllReadInform << SAMFileNumber << '\t' << Read.first;
                            for (int j = 0; j < Read.second.size(); j++) {
                                if (j == 1) {
                                    AllReadInform << '\t' << Group_index;
                                } else if (j == 2) {
                                    AllReadInform << '\t' << Group_low;
                                } else if (j == 3) {
                                    AllReadInform << '\t' << Group_high;
                                } else {
                                    AllReadInform << '\t' << Read.second[j];
                                }
                            }
                            chr_name = Read.second[0];
                            AllReadInform << '\n';
                        }

                        Readpos = AllReadInform.tellp();
                        Chunk_Bang.coverage2pos.insert({{Group_low, Group_high}, {Chunk_Bang.reads_pointer.back(),Readpos}});
                        Chunk_Bang.reads_pointer.push_back(Readpos); 
                        Chunk_Bang.reads_pointer.push_back(Readpos);  
                        Chunk_Bang.group_reads_number.push_back(Group_temporary.size());                   
                        
                        auto it = Chunk_Bang.chr_coverage.find(chr_name);
                        if (it != Chunk_Bang.chr_coverage.end()){
                            Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                        } else {
                            Chunk_Bang.chr_coverage[chr_name] = {};
                            Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                        }
                   
                        Group_temporary.clear();
                        Group_temporary[This_Line.read_name] = This_Line.tokens;    
                        Group_high = this_coverage[1];
                        Group_low = this_coverage[0];                                             
                    }
                } else {
                    Group_temporary[This_Line.read_name] = This_Line.tokens;
                    Group_high = this_coverage[1];
                    Group_low = this_coverage[0];
                }
            }
            earlyGroup = thisGroup; 
            last_chr = This_Line.tokens[0];
        }

        if (file_number == (Read_x_vec.size()-1) && SmallSamFile.eof()){
            Group_index++;
            for (const auto& Read:Group_temporary){
                AllReadInform << SAMFileNumber << '\t' << Read.first;
                for (int j = 0; j < Read.second.size(); j++) {
                    if (j == 1) {
                        AllReadInform << '\t' << Group_index;
                    } else if (j == 2) {
                        AllReadInform << '\t' << Group_low;
                    } else if (j == 3) {
                        AllReadInform << '\t' << Group_high;
                    } else {
                        AllReadInform << '\t' << Read.second[j];
                    }
                }
                chr_name = Read.second[0];
                AllReadInform << '\n';
            }
            Readpos = AllReadInform.tellp();
            Chunk_Bang.coverage2pos.insert({{Group_low, Group_high}, {Chunk_Bang.reads_pointer.back(),Readpos}});
            Chunk_Bang.reads_pointer.push_back(Readpos);
            Chunk_Bang.group_reads_number.push_back(Group_temporary.size());

            auto it = Chunk_Bang.chr_coverage.find(chr_name);
            if (it != Chunk_Bang.chr_coverage.end()) {
                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
            } else {
                Chunk_Bang.chr_coverage[chr_name] = {};
                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
            }
        }
        SmallSamFile.close();
        std::remove(fileName.c_str());
    }
    AllReadInform.close();
    Chunk_Bang.FileNo = 1;
    return Chunk_Bang;
}


std::map<std::string, std::vector<std::array<int,2>>> Merge_Read_Interval(std::map<int, std::map<std::string, std::vector<std::array<int,2>>>>& FileCoverage) {

    std::map<std::string, std::vector<std::array<int,2>>> chr_range_temp;
    std::vector<std::array<int,2>> raw_intervals;
    std::vector<std::array<int,2>> sorted_intervals;
    if (FileCoverage.size() != 0) {

        for (const auto& eachFile:FileCoverage) {
            if (eachFile.second.size() != 0) {

                for (const auto& eachChr:eachFile.second) {
                    auto itchr = chr_range_temp.find(eachChr.first);
                    if (itchr != chr_range_temp.end()) {
                        chr_range_temp[eachChr.first].insert(chr_range_temp[eachChr.first].end(), eachChr.second.begin(), eachChr.second.end());
                    } else {
                        chr_range_temp[eachChr.first] = {};
                        chr_range_temp[eachChr.first].insert(chr_range_temp[eachChr.first].end(), eachChr.second.begin(), eachChr.second.end());
                    }
                }
            }
        }
    }

    std::map<std::string, std::vector<std::array<int,2>>> ChrInterval;
    if (chr_range_temp.size() != 0) {
        for (const auto& eachChr:chr_range_temp) {
            raw_intervals = eachChr.second;
            sort(raw_intervals.begin(), raw_intervals.end());
            sorted_intervals.clear();

            for (int i=0; i<raw_intervals.size();) {
                int left = raw_intervals[i][0];
                int right = raw_intervals[i][1];
                while (i+1 < raw_intervals.size() && right >= raw_intervals[i+1][0]){
                    right = std::max(right, raw_intervals[i+1][1]);
                    i++;
                }
                sorted_intervals.push_back({left, right});
                i++;
            }
            ChrInterval[eachChr.first] = sorted_intervals;
        }
    }
    return ChrInterval;
}




void Merge_Read_Big_Files_Group(const int& NEWNEW,
                                const std::string& output_path, 
                                const std::array<int,2>& Regin_BE,
                                std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>>& file_pointers,
                                const std::string& chrchr,
                                std::vector<std::streampos>& group_read_pointer,
                                std::vector<int>& group_read_number,
                                std::ofstream& BigFile) {

    std::streampos ThisFile_ThisGroup_pointer;
    std::streampos end;
    std::string line;
    std::map<int, std::vector<std::vector<std::string>>> Group_store;

    std::map<int, std::array<std::streampos,2>> This_Coverage_pointers = file_pointers[chrchr][Regin_BE];

    if (This_Coverage_pointers.size() > 0) {
        for (auto it = This_Coverage_pointers.begin(); it != This_Coverage_pointers.end(); ++it) {
            int i = it->first;
            Group_store[i] = {};
            const std::array<std::streampos, 2>& positions = it->second;

            std::string ReadFilePath = output_path + "/sam_" + std::to_string(i) + "/All_Read.txt";

            std::ifstream samfile_one(ReadFilePath);

            ThisFile_ThisGroup_pointer = positions[0];
            end = positions[1];
            samfile_one.seekg(ThisFile_ThisGroup_pointer);
            while (getline(samfile_one, line)) {
                Split_Result This_Line = get_sj_split(line, '\t');
                Group_store[i].push_back(This_Line.tokens);
                ThisFile_ThisGroup_pointer = samfile_one.tellg();
                if (ThisFile_ThisGroup_pointer < end) {
                    continue;
                } else {
                    break;
                }
            }
            samfile_one.close();
        }
    }

    int count = 0;
    int line_count = -1;
    if (BigFile.is_open()) 
    {
        {
            std::unique_lock<std::mutex> lock(bigMutex);
            for (const auto& eachFile:Group_store) {
                for (const auto& eachLine:eachFile.second) {
                    count++;
                    line_count = -1;
                    for (const auto& eachpart:eachLine) {
                        line_count = line_count + 1;
                        if (line_count == 0) {
                            BigFile << eachpart;
                            continue;
                        } else if (line_count == 3) {
                            BigFile << '\t' << NEWNEW;
                            continue;
                        } else if (line_count == 4) {
                            BigFile << '\t' << Regin_BE[0];
                            continue;
                        } else if (line_count == 5) {
                            BigFile << '\t' << Regin_BE[1];
                            continue;
                        } else {
                            BigFile << '\t' << eachpart;
                        }
                    }
                    BigFile << '\n';
                }
            }
            group_read_number.push_back(count);
            ThisFile_ThisGroup_pointer = BigFile.tellp();
            group_read_pointer.push_back(ThisFile_ThisGroup_pointer); 
            group_read_pointer.push_back(ThisFile_ThisGroup_pointer); 
        }
    } else {
        std::cout << "^-^ *** Merge large files failed *** ^-^" << std::endl;
    }
}


std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> get_pointers(
        std::map<std::string, std::vector<std::array<int,2>>>& chrcoverage,
        std::map<int, std::map<std::string, std::vector<std::array<int,2>>>>& FileCoverage,
        std::map<int, std::map<std::array<int,2>, std::array<std::streampos,2>>>& groups_pointer) {

    std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> chr_pointer_coverage;

    for (const auto& EachChr:chrcoverage) {
        chr_pointer_coverage[EachChr.first] = {};
        for (const auto& Region:EachChr.second) {
            chr_pointer_coverage[EachChr.first][Region] = {};
        }
    }

    std::array<std::streampos,2> thisFile_group_be;

    for (const auto& EachFile:FileCoverage) {

        int FileNumber = EachFile.first;
        std::map<std::array<int,2>, std::array<std::streampos,2>> thisFile_group_pointer = groups_pointer[FileNumber];

        for (const auto& EachChr:EachFile.second) {
            int Rhino = 0; 
            int Rhino_small = 0;
            std::vector<std::streampos> temp_pointer; 
            std::string ChrName = EachChr.first;
            std::vector<std::array<int,2>> Big_Rhino = chrcoverage[ChrName]; 

            for (const auto& small:EachChr.second) {

                Rhino_small = Rhino_small + 1;

                for (int i = Rhino; i < Big_Rhino.size(); i++) {
                    if (small[0] >= Big_Rhino[i][0] && small[1] <= Big_Rhino[i][1]) {

                        temp_pointer.push_back(thisFile_group_pointer[small][0]);
                        temp_pointer.push_back(thisFile_group_pointer[small][1]);
 
                        if (Rhino_small == EachChr.second.size()) {
                            thisFile_group_be[0] = temp_pointer[0];
                            thisFile_group_be[1] = temp_pointer[temp_pointer.size()-1];
                            chr_pointer_coverage[ChrName][Big_Rhino[i]].insert({FileNumber, thisFile_group_be});
                            temp_pointer.clear();
                        }
                        Rhino = i;
                        break;
                    } else if (small[0] > Big_Rhino[i][1]) {

                        if (temp_pointer.size() > 0) {
                            thisFile_group_be[0] = temp_pointer[0];
                            thisFile_group_be[1] = temp_pointer[temp_pointer.size()-1];
                            temp_pointer.clear();
                            chr_pointer_coverage[ChrName][Big_Rhino[i]].insert({FileNumber, thisFile_group_be});
                        }
                    } else {

                        continue;
                    }

                }
            }
        }

    }
    return chr_pointer_coverage;
}


FileSplit thread_all_read_sam_files(const std::string& sam_file_path, 
                                    std::vector<std::string>& samfilevec,
                                    const int& numThreads, 
                                    const std::string& outputPath, 
                                    const int& SJDistance, 
                                    std::unordered_map<std::string, std::string>& fasta_ref) {
    std::vector<std::string> sam_file_vec;
    if (samfilevec.size() == 1) {
        sam_file_vec.push_back(sam_file_path);
    } else {
        for (const auto& eachFile:samfilevec) {
            std::string chunkFilePath = joinPath(sam_file_path, eachFile);
            sam_file_vec.push_back(chunkFilePath);
        }
    }

    FileSplit BigBang;
    std::map<int, std::map<std::string, std::vector<std::array<int,2>>>> File_chr_coverage;
    std::map<int, std::map<std::array<int,2>, std::array<std::streampos,2>>> File_group_pointer;
    std::vector<std::streampos> startposVec;
    std::vector<std::streampos> endposVec;
    if (sam_file_vec.size() != 0) {

        ThreadPool Bigfilepool(numThreads);
        std::vector<std::future<void>> myJobs;
        for (int samFileNumber = 0; samFileNumber < sam_file_vec.size(); ++samFileNumber) {
            
            startposVec.clear();
            endposVec.clear();

            std::cout << "******* " << "Start processing SAM File " << samFileNumber << " *******" << std::endl;
            std::string chunkFilePath = "sam_" + std::to_string(samFileNumber);
            chunkFilePath = joinPath(outputPath, chunkFilePath);

            if (!directoryExists(chunkFilePath)) {
                createDirectory(chunkFilePath);      
            }
            struct stat statBuf;
            if (stat(sam_file_vec[samFileNumber].c_str(), &statBuf) != 0){
                std::cout << "We can't get file size !" << sam_file_vec[samFileNumber] << std::endl;
                exit(EXIT_FAILURE);
            }
            std::streampos fileSize = statBuf.st_size;

            std::streampos chunkSize = fileSize / numThreads;

            std::ifstream one_sam(sam_file_vec[samFileNumber]);

            if (!one_sam.is_open())	{
                std::cout << "We can't open SAM file !" << sam_file_vec[samFileNumber] << std::endl;
                exit(EXIT_FAILURE);
            }

            for (int i = 0; i < numThreads; ++i) {
                std::streampos startPos = i * chunkSize;
                std::streampos endPos = (i == numThreads - 1) ? (fileSize) : static_cast<std::streampos>((i + 1) * chunkSize);

                if (i > 0) {
                    startPos = findNextLineStart(one_sam, startPos);
                }
                startposVec.push_back(startPos);
                endposVec.push_back(endPos);
            }

            one_sam.close();
  
            myJobs.clear();
            for (int i = 0; i < numThreads; ++i) {
                myJobs.emplace_back(Bigfilepool.enqueue([&,i]() {
                    processChunk(
                        sam_file_vec[samFileNumber],
                        startposVec[i], 
                        endposVec[i],
                        chunkFilePath, 
                        i, 
                        SJDistance, 
                        std::ref(fasta_ref)
                    );
                }));
            }
            for (auto& future : myJobs) {
                future.get();  // 等待每个任务完成
            }
            std::cout << "^-^ [" << sam_file_vec[samFileNumber] << "] All threads are finished generating small files! ^-^" << std::endl;
            
            std::cout << "^-^ Start of merge small files ! ^-^" << std::endl;
            BigBang = Merge_Read_Small_Files(chunkFilePath, samFileNumber);
            std::cout << "^-^ End of merge small files ! ^-^" << std::endl;
            
            if (sam_file_vec.size() > 1) {

                File_chr_coverage[samFileNumber] = BigBang.chr_coverage;
                File_group_pointer[samFileNumber] = BigBang.coverage2pos;

            }

        }
        if (sam_file_vec.size() > 1) {
            std::cout << "^-^ The number of sam files is greater than one. Large files need to be merged. ^-^" << std::endl;

            std::map<std::string, std::vector<std::array<int,2>>> ChrCoverage = Merge_Read_Interval(File_chr_coverage);


            std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> ChrCoverage_SmallPointer = 
                                                                        get_pointers(ChrCoverage, File_chr_coverage, File_group_pointer);
   
            std::vector<std::streampos> Allreads_pointer;
            std::vector<int> Allreads_group_number;            
            // 输出的文件名;
            std::string FinallyFile = joinPath(outputPath, "MergedRead.txt");
            std::ofstream FinallyReadInform(FinallyFile, std::ios::trunc); 
            std::streampos Readpos = FinallyReadInform.tellp();
            Allreads_pointer.push_back(Readpos);


            myJobs.clear();
            int new_group = 0;

            std::cout << "^-^ Start of merge large File ! ^-^" << std::endl;

            for (const auto& eachChr:ChrCoverage_SmallPointer) {
                std::string ChrName = eachChr.first;
                // std::cout << ChrName << std::endl;
                for (const auto& eachRegion:eachChr.second) {
                    std::array<int,2> ourRegion = eachRegion.first;
                    new_group = new_group + 1;
                    myJobs.emplace_back(Bigfilepool.enqueue([new_group, ourRegion, ChrName, &outputPath, &ChrCoverage_SmallPointer, &Allreads_pointer, &Allreads_group_number, &FinallyReadInform]() {
                        Merge_Read_Big_Files_Group(
                            new_group,
                            outputPath,
                            ourRegion,
                            ChrCoverage_SmallPointer,
                            ChrName,
                            Allreads_pointer,
                            Allreads_group_number,
                            FinallyReadInform
                        );
                    }));
                }
            }
            for (auto& future : myJobs) {
                future.get();
            }

            std::cout << "^-^ End of merge large File ! ^-^" << std::endl;

            BigBang.readtxt_path = FinallyFile;
            BigBang.group_reads_number = Allreads_group_number;
            BigBang.reads_pointer = Allreads_pointer;
            BigBang.FileNo = sam_file_vec.size();

            FinallyReadInform.close();
        }
    } else {
        std::cerr << "* There are no files in the folder! *" << strerror(errno) << std::endl;
    }
    return BigBang;
}


struct Line_Split{
    std::string read_file; //0
    std::string read_name; //1
    std::string chr_name; //2
    std::string Group_index; //3
    std::array<int,2> Group_coverage; //4,5
    std::array<int,2> read_coverage; //6,7
    int read_length; // 这个就是追踪文件的时候可以画Coverage; //8
    std::vector<std::array<int,2>> read_SJ;
    std::vector<int> read_sj_quality;
    
};

Line_Split MakeLineSplit(const std::string& s, char delimiter){
    // std::cout << s << std::endl;
    // 初始化返回值结构体;
    Line_Split read_result = {};
    // 字符串流;
    std::istringstream iss(s);
    // 定义字符串流中的每一段字符串;
    std::string token;
    // std::cout << iss << std::endl;
    int number = 0;
    // 用于二次处理;
    std::vector<int> secondProcess;
    std::array<int,2> Sj;
    while (std::getline(iss, token, delimiter)){
        number = number + 1;
        if (number == 1) {
            read_result.read_file = token; 
        } else if (number == 2) {
            read_result.read_name = token;
        } else if (number == 3) {
            read_result.chr_name = token;
        } else if (number == 4) {
            read_result.Group_index = token;
        } else if (number == 5) {
            read_result.Group_coverage[0] = std::stoi(token);
        } else if (number == 6) {
            read_result.Group_coverage[1] = std::stoi(token);
        } else if (number == 7) {
            read_result.read_coverage[0] = std::stoi(token);
        } else if (number == 8) {
            read_result.read_coverage[1] = std::stoi(token);
        } else if (number == 9) {
            read_result.read_length = std::stoi(token);
        } else if (number > 9) {
            secondProcess.push_back(std::stoi(token));
        }
    }
    // std::cout << std::endl;
    // 两次处理;
    for (size_t i = 0; i < secondProcess.size(); i += 3) {
        Sj = {secondProcess[i], secondProcess[i+1]};
        read_result.read_SJ.push_back(Sj);
        read_result.read_sj_quality.push_back(secondProcess[i+2]);
    }
    return read_result;
}



struct GroupInformation {
    std::string GroupIndex;
    std::string chrName; //染色体;
    std::array<int,2> GroupCoverage; //Group范围;
    std::unordered_map<std::string, std::string> GroupReadFiles; //Group每条reads的所属文件; (公共, single exon 与 Multi exon都有);
    std::unordered_map<std::string, std::vector<std::array<int,2>>> GroupReadSjs; //Group每条reads的SJs;
    std::unordered_map<std::string, std::array<int,2>> GroupReadCoverage; //Group每条reads的开头结尾;
    std::map<std::array<int,2>, int> GroupSjs; //Group中根据Reads得到的高置信SJ;
    std::unordered_map<std::string, std::vector<int>> GroupSigns; //Group中所有reads的信号;
    std::unordered_map<std::string, std::array<int,2>> GroupSingleExon; // single exon的reads;
};

GroupInformation knowGroupInformation(std::streampos& startpos, std::streampos& endpos, const std::string& sam_file_path, const int& Sj_Support_Read_Number) {
    GroupInformation groupinformation;
    std::ifstream FinalSamFile(sam_file_path);
    FinalSamFile.seekg(startpos);
    // 记录当前, 也就是文件刚开始的位置;
    std::streampos Current_Position = FinalSamFile.tellg();
    std::streampos Last_Position = Current_Position;
    std::string line;
    std::map<std::array<int,2>, std::array<int,2>> Temp_Sjs;
    // 读取整个Group;
    while (getline(FinalSamFile, line)) {
        // 记录处理到的位置;
        Last_Position = Current_Position;
        Current_Position = FinalSamFile.tellg();
        if (Last_Position < endpos || Last_Position == FinalSamFile.eof()) {
            // 拆解每一行;
            Line_Split Line_result = MakeLineSplit(line, '\t');
            groupinformation.GroupIndex = Line_result.Group_index;
            groupinformation.chrName = Line_result.chr_name;
            groupinformation.GroupCoverage = Line_result.Group_coverage;
            groupinformation.GroupReadFiles[Line_result.read_name] = Line_result.read_file;
            groupinformation.GroupReadCoverage[Line_result.read_name] = Line_result.read_coverage;
            if (Line_result.read_SJ.size() > 0) {
                groupinformation.GroupReadSjs[Line_result.read_name] = Line_result.read_SJ;
                groupinformation.GroupSigns[Line_result.read_name] = Line_result.read_sj_quality;
                // 将所有的Sj循环遍历填充到Temp_Sjs中;
                for (int i = 0; i < Line_result.read_sj_quality.size(); i++) {
                    if (Line_result.read_sj_quality[i] > 0) {
                        std::array<int,2> Temp_array = Line_result.read_SJ[i];
                        auto it = Temp_Sjs.find(Temp_array);
                        if (it != Temp_Sjs.end()) {
                            Temp_Sjs[Temp_array][0] = Temp_Sjs[Temp_array][0] + 1;
                        } else {
                            Temp_Sjs[Temp_array] = {1, Line_result.read_sj_quality[i]};
                        }
                    }
                }
            } else {
                groupinformation.GroupSingleExon[Line_result.read_name] = Line_result.read_coverage;
            }
        } else {
            break;
        }
    }
    FinalSamFile.close();
    // 整个cluster结束之后, 撰写结果;
    for (const auto& Sj:Temp_Sjs) {
        if (Sj.second[0] > Sj_Support_Read_Number) {
            groupinformation.GroupSjs[Sj.first] = Sj.second[1];
        }
    }
    return groupinformation;
}



struct GroupAnnotation
{
    std::unordered_map<std::string, std::vector<std::array<int,2>>> Group_Annotations;
    std::set<std::string> Group_GeneSet;
};

GroupAnnotation get_group_annotation(std::unordered_map<std::string, std::array<int,2>>& AnnoCoverage, std::unordered_map<std::string, std::vector<std::array<int,2>>>& AnnoisoSJ, std::array<int,2>& Group_Coverage){
    //输出的数据结构就是unordered_map;
    GroupAnnotation ThisGroupAnnotations;
    std::string extracted;

    // int a = 0;
    //对当前的注释遍历;
    for (const auto& EachAnno:AnnoCoverage) {
        if ((EachAnno.second[0] >= Group_Coverage[0]) && (EachAnno.second[1] <= Group_Coverage[1])){
            //满足条件的存储一下;
            ThisGroupAnnotations.Group_Annotations[EachAnno.first] = AnnoisoSJ[EachAnno.first];
            // a++;
            size_t pos = EachAnno.first.find("|");
            if (pos != std::string::npos) {
                extracted = EachAnno.first.substr(0, pos);
            }
            ThisGroupAnnotations.Group_GeneSet.insert(extracted);
        }
    }
    return ThisGroupAnnotations;
}

// 得到Group范围内的单个exon的注释;
std::unordered_map<std::string, std::array<int,2>> get_group_single_exon_annotation(
                                    std::unordered_map<std::string, std::array<int,2>>& AnnoSE, 
                                    std::array<int,2>& Group_Coverage) {
    //输出的数据结构就是unordered_map;
    std::unordered_map<std::string, std::array<int,2>> ThisGroupSingleExonAnnotations;
    //对当前的注释遍历;
    for (const auto& EachAnno:AnnoSE) {
        if ((EachAnno.second[0] >= Group_Coverage[0]) && (EachAnno.second[1] <= Group_Coverage[1])){
            //满足条件的存储一下;
            ThisGroupSingleExonAnnotations[EachAnno.first] = EachAnno.second; 
        }
    }
    return ThisGroupSingleExonAnnotations;
}

// 得到Group范围内的单个exon及其对应的reads;
std::unordered_map<std::string, std::vector<std::string>> get_group_single_exon_reads(
                                    std::unordered_map<std::string, std::array<int,2>>& ThisGroupSingleExonAnnotations,
                                    std::unordered_map<std::string, std::array<int,2>>& GroupSingleExonReads,
                                    int& Edge) {
    // 输出这个注释对应的reads;
    std::unordered_map<std::string, std::vector<std::string>> SingleReads;
    // 初始化;
    for (const auto& eachAnno:ThisGroupSingleExonAnnotations) {
        SingleReads[eachAnno.first] = {};
    }
    // 开始对每个single exon 循环;
    for (const auto& EachSing:GroupSingleExonReads) {
        // 每个注释;
        for (const auto& EachAnno:ThisGroupSingleExonAnnotations) {
            if ((EachSing.second[1]-EachAnno.second[1]<Edge) && (EachAnno.second[0]-EachSing.second[0]<Edge)) {
                SingleReads[EachAnno.first].push_back(EachSing.first);
                break;
            }
        }
    }
    return SingleReads;
}


std::unordered_map<std::string, std::map<std::string, double>> write_single_exon_gtf_trace(int& FileNo,
    std::unordered_map<std::string, std::array<int,2>>& groupsingleexon,
    std::unordered_map<std::string, std::string>& groupreadfiles,
    std::unordered_map<std::string, std::vector<std::string>>& singleexonwithreads,
    std::unordered_map<std::string, std::array<int,2>>& singleexongroupannotation,
    std::ofstream& Updated_Files, std::ofstream& Trace, std::ofstream& isoformFilePath,
    std::unordered_map<std::string, std::string>& GTF_Transcript_Strand,
    std::string& chrname) {
    // 最终的结果还是return出去吧;
    // 最后一起写出来;
    // 这个奇妙了, 我可以unmap<isoform:<file:counts>>
    std::unordered_map<std::string, std::map<std::string, double>> FileSingleExonNumber;
    std::unordered_map<std::string, std::map<std::string, double>> FileSingleGeneNumber;
    std::string first_part;
    std::string second_part;
    // 以上初始化最终结果;
    // 首先看看有几个文件哈;
    if (FileNo > 1) {
        if (Updated_Files.is_open()) {
            // 多个sam文件;
            std::map<std::string, double> thisIsoform_FileCounts;
            // 首先遍历已经根据注释分好的组;
            if (singleexonwithreads.size() > 0) {
                for (const auto& eachAnno:singleexonwithreads) {
                    // 可以写每个转录本的gtf;
                    if (eachAnno.second.size() > 0) {
                        size_t pos = eachAnno.first.find('|');
                        if (pos != std::string::npos) {
                            // 根据 "|" 将字符串分割为两部分
                            first_part = eachAnno.first.substr(0, pos);
                            second_part = eachAnno.first.substr(pos + 1);
                        }
                        {
                            std::unique_lock<std::mutex> lock(updatedGtfMutex);
                            Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t' << singleexongroupannotation[eachAnno.first][0] << '\t' << singleexongroupannotation[eachAnno.first][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[eachAnno.first] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\";" << '\n';
                            Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t' << singleexongroupannotation[eachAnno.first][0] << '\t' << singleexongroupannotation[eachAnno.first][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[eachAnno.first] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\"; exon_number \"" << 0 << "\";" << '\n';
                        }
                        // 每个注释转录本开始时的归零;
                        thisIsoform_FileCounts.clear();
                        for (int i = 0; i < FileNo; i++) {
                            thisIsoform_FileCounts[std::to_string(i)] = 0;
                        }
                        // 开始数数;
                        for (const auto& eachRead:eachAnno.second) {
                            std::string thisReadFile = groupreadfiles[eachRead];
                            thisIsoform_FileCounts[thisReadFile] = thisIsoform_FileCounts[thisReadFile] + 1;
                            // 下面可以写每个read的Trace;
                            {
                                std::unique_lock<std::mutex> lock(traceMutex);
                                Trace << eachRead << '\t' << "single_exon" << '\t' << second_part << '\t' << first_part << '\t' << thisReadFile << '\n'; 
                            }
                        }
                        // 最终数量结果写入;
                        FileSingleExonNumber[eachAnno.first] = thisIsoform_FileCounts;
                        // 一个isoform完毕;
                    }
                }
            }            
        }
    } else {
        if (Updated_Files.is_open()) {
            // 一个sam文件;
            // 首先遍历已经根据注释分好的组;
            if (singleexonwithreads.size() > 0) {
                for (const auto& eachAnno:singleexonwithreads) {
                    // 可以写每个转录本的gtf;
                    if (eachAnno.second.size() > 0) {
                        size_t pos = eachAnno.first.find('|');
                        if (pos != std::string::npos) {
                            // 根据 "|" 将字符串分割为两部分
                            first_part = eachAnno.first.substr(0, pos);
                            second_part = eachAnno.first.substr(pos + 1);
                        }
                        {
                            std::unique_lock<std::mutex> lock(updatedGtfMutex);
                            Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t' << singleexongroupannotation[eachAnno.first][0] << '\t' << singleexongroupannotation[eachAnno.first][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[eachAnno.first] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\";" << '\n';
                            Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t' << singleexongroupannotation[eachAnno.first][0] << '\t' << singleexongroupannotation[eachAnno.first][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[eachAnno.first] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\"; exon_number \"" << 0 << "\";" << '\n';
                        }
                        // 开始数数;
                        FileSingleExonNumber[eachAnno.first]["0"] = eachAnno.second.size();
                        for (const auto& eachRead:eachAnno.second) {
                            // 下面可以写每个read的Trace;
                            {
                                std::unique_lock<std::mutex> lock(traceMutex);
                                Trace << eachRead << '\t' << "single_exon" << '\t' << second_part << '\t' << first_part << '\t' << 0 << '\n'; 
                            }
                        }
                        // 一个isoform完毕;
                    }

                }
            }                        
        }
    }
    // 输出单个exon的定量, 并计算出基因的数量;
    for (const auto& eachAnno:FileSingleExonNumber) {
        size_t pos = eachAnno.first.find('|');
        if (pos != std::string::npos) {
            // 根据 "|" 将字符串分割为两部分
            first_part = eachAnno.first.substr(0, pos);
            second_part = eachAnno.first.substr(pos + 1);
        }
        // 写到文件里;
        {
            std::unique_lock<std::mutex> lock(isoformCountMutex);
            isoformFilePath << second_part << '\t' << first_part;
            for (const auto& eachFile:eachAnno.second) {
                isoformFilePath << '\t' << eachFile.second;
            }
            isoformFilePath << '\n';
        }
        // 生成基因的变量;
        auto it = FileSingleGeneNumber.find(first_part); //基因;
        if (it != FileSingleGeneNumber.end()) {
            // 已经有了这个基因;
            // 叠加就可以;
            for (const auto& eachFile:eachAnno.second) {
                FileSingleGeneNumber[first_part][eachFile.first] = FileSingleGeneNumber[first_part][eachFile.first] + eachFile.second;
            }
        } else {
            // 尚且没有这个基因;
            FileSingleGeneNumber[first_part] = {};
            for (const auto& eachFile:eachAnno.second) {
                FileSingleGeneNumber[first_part][eachFile.first] = eachFile.second;
            }
        }
    }
    // return FileSingleExonNumber;
    return FileSingleGeneNumber;
}



void merge_high_sj(std::unordered_map<std::string, std::vector<std::array<int,2>>>& group_annotations, std::unordered_map<std::string, std::string>& transcript_strand, std::map<std::array<int,2>, int>& Grouphighsjs) {
    // 注释的所有sj, 所在的链, 高置信区间;
    if (group_annotations.size() != 0) {
        for (const auto& isoform:group_annotations) {
            for (const auto& sj:isoform.second) {
                auto it = Grouphighsjs.find(sj);
                if (it == Grouphighsjs.end()) {
                    if (transcript_strand[isoform.first] == "+") {
                        Grouphighsjs[sj] = 1;
                    } else {
                        Grouphighsjs[sj] = 2;
                    }
                }
            }
        }
    } else {
        return;
    }
}


std::unordered_map<std::size_t, std::vector<std::string>> classifyReadsVec(std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs) {
    std::unordered_map<std::size_t, std::vector<std::string>> ClustersReads; //用于存储cluster;
    std::map<std::vector<std::array<int,2>>, std::size_t> readsj2sizet;
    //提取每个map中存储的每个reads的SJ的小区间;
    std::vector<std::array<int,2>> vec;
    std::string vecname;
    std::size_t cluster = 0;
    std::size_t clusterNumber = 0;
    for (const auto& pair : groupreadsjs) {
        vecname = pair.first;
        vec = pair.second;
        auto itsj = readsj2sizet.find(vec);
        if (itsj != readsj2sizet.end()) {
            cluster = readsj2sizet[vec];
            ClustersReads[cluster].push_back(vecname);
        } else {
            // 新的cluster;
            readsj2sizet[vec] = clusterNumber;
            ClustersReads[clusterNumber].push_back(vecname);
            clusterNumber = clusterNumber + 1;
        }
    }
    return ClustersReads;
}


std::unordered_map<std::size_t, std::array<int,2>> get_every_cluster_begin_end(std::unordered_map<std::size_t, std::vector<std::string>>& ClusterReads, 
                                                                               std::unordered_map<std::string, std::array<int,2>>& AllReadBE){
    std::unordered_map<std::size_t, std::array<int,2>> CluserBeginEnd;
    std::vector<int> ExonbeginVec;
    std::vector<int> ExonendVec;
    for (const auto& Tuan:ClusterReads){
        ExonbeginVec.clear();
        ExonendVec.clear();
        CluserBeginEnd[Tuan.first] = {};
        for (const auto& oneread:Tuan.second){
            ExonbeginVec.push_back(AllReadBE[oneread][0]);
            ExonendVec.push_back(AllReadBE[oneread][1]);
        }
        //分别对两个变量取值;
        auto min_it = std::min_element(ExonbeginVec.begin(), ExonbeginVec.end());
        auto max_it = std::max_element(ExonendVec.begin(), ExonendVec.end());
        CluserBeginEnd[Tuan.first][0] = *min_it;
        CluserBeginEnd[Tuan.first][1] = *max_it;
    }
    return CluserBeginEnd;
}


struct ReferenceCluster
{
    std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>> FSM;
    std::map<std::size_t, std::vector<std::string>> ISM; //保存FSM的ISM和Others的ISM;
    std::unordered_map<std::size_t, std::vector<std::string>> Others;
};

ReferenceCluster get_FSM_and_others(std::unordered_map<std::size_t, std::vector<std::string>>& ClustersReads, 
                                   std::unordered_map<std::string, std::vector<std::array<int,2>>>& GroupAnno, 
                                   std::unordered_map<std::string, std::array<int,2>>& GroupAnnoBE, 
                                   std::unordered_map<std::string, std::array<int,2>>& AllReadBE, 
                                   std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs,
                                   std::unordered_map<std::string, std::string>& GroupReadFile, 
                                   std::ofstream& Trace){
    // std::mutex& TraceMutex;
    ReferenceCluster FSMISMO = {};
    //首先对所有的高质量的cluster遍历;
    //SIRV可以看所在的是哪个基因;
    std::vector<std::array<int,2>> each_read_SJs; //reads的所有SJ;

    //新生成的变量;
    std::vector<std::string> statis_FSM; //满足条件的FSM的名称;
    std::vector<std::string> statis_ISM; //满足条件的ISM的名称;
    std::vector<int> ISM_flag; //ISM中每个SJ的索引;
    // int count = 0; //满足条件的ISM个数;
    std::string FSM_name; //保存输出变量时FSM的名称;
    std::vector<int> minFSM; //每个FSM与reads的距离;
    int minmin = 0; //FSM最小位置的索引;

    std::string first_part;
    std::string second_part; 

    //遍历cluster;
    for (auto it = ClustersReads.begin(); it != ClustersReads.end(); ++it) {
        //每一项的key表示这个cluster的SJ; value表示所有的reads的名称;
        each_read_SJs = AllSJs[(it->second)[0]];
        
        //遍历注释;
        statis_FSM.clear(); //判断FSM的, 判断多少个, 一般情况下就是1个;
        statis_ISM.clear();
        for (const auto& pair : GroupAnno){
            //pair.first表示基因名称如SIRV101, pair.second表示注释的所有SJ;
            if (pair.second == each_read_SJs){
                //这个判断条件表示 注释的SJ是否等于cluster的SJ;
                statis_FSM.push_back(pair.first);
            } else {
                //不相等的情况;
                //注释SJ的个数 必须大于 clusterSJ的个数;
                if (pair.second.size() > each_read_SJs.size()){
                    //接下来就看看是不是ISM; 不需要知道是谁的ISM, 后续还会追究, 所以通过数字就可以了解;
                    ISM_flag.clear();
                    
                    for (const auto& element : each_read_SJs) {
                        int it_index = std::find(pair.second.begin(), pair.second.end(), element) - pair.second.begin();
                        if (it_index == pair.second.size()){
                            ISM_flag.push_back(10000);   
                            break;
                        } else{
                            ISM_flag.push_back(it_index);
                        }
                    }
                    if ((ISM_flag[ISM_flag.size()-1] != 10000) && (ISM_flag[0] != 10000) && (ISM_flag[ISM_flag.size()-1] - ISM_flag[0] == ISM_flag.size() - 1)) {statis_ISM.push_back(pair.first);}
                }
                //这里说明是这个cluster, 比当前的注释转录本的SJ个数要多;
            }
        }
        //遍历一遍注释后;继续看满足条件的个数;
        //先看FSM;
        if (statis_FSM.size() == 1){
            FSM_name = statis_FSM[0];
            FSMISMO.FSM[FSM_name].first = each_read_SJs;
            FSMISMO.FSM[FSM_name].second.insert(FSMISMO.FSM[FSM_name].second.end(), ClustersReads[it->first].begin(), ClustersReads[it->first].end());
            //将这个cluster中的每一条reads都输出到追踪文件;
            size_t pos = FSM_name.find('|');
            if (pos != std::string::npos){
                first_part = FSM_name.substr(0, pos);
                second_part = FSM_name.substr(pos + 1);
            }
            if (Trace.is_open()){
                for (const auto& EachRead:ClustersReads[it->first]){
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << EachRead << '\t' << "FSM" << '\t' << second_part << '\t' << first_part << '\t' << GroupReadFile[EachRead] << '\n'; 
                }
            }                        

        } else if (statis_FSM.size() > 1){
            
            //这种情况就需要对每条reads分别看开头和结尾;
            //初始化将数量标记为0;
            for (const auto& refx : statis_FSM){
                FSM_name = refx;
                FSMISMO.FSM[FSM_name].first = each_read_SJs;
                FSMISMO.FSM[FSM_name].second = {};
            }
            for (const auto& readx : ClustersReads[it->first]) { //readx是每条read的名称;
                minFSM.clear();
                for (const auto& refx : statis_FSM){
                    minFSM.push_back(std::abs(AllReadBE[readx][0] - GroupAnnoBE[refx][0]) + std::abs(AllReadBE[readx][1] - GroupAnnoBE[refx][1]));
                }
                //每条reads的所有候选FSM计算完毕;
                minmin = std::min_element(minFSM.begin(), minFSM.end()) - minFSM.begin();
                //因此这条reads的归于索引为minmin的isoform;
                FSM_name = statis_FSM[minmin];
                FSMISMO.FSM[FSM_name].second.push_back(readx); //有个小bug, 如果都一样的话, 只会给前一个;
                //将这个cluster中的每一条reads都输出到追踪文件;
                size_t pos = FSM_name.find('|');
                if (pos != std::string::npos){
                    first_part = FSM_name.substr(0, pos);
                    second_part = FSM_name.substr(pos + 1);
                }
                {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << readx << '\t' << "FSM" << '\t' << second_part << '\t' << first_part << '\t' << GroupReadFile[readx] << '\n'; 
                }                
            } //分配结束;

        } else{ //只有FSM没有的情况下, 才继续看ISM;
            //然后看ISM;
            if (statis_ISM.size() == 0){
                //这种的情况是其他类;
                FSMISMO.Others[it->first] = ClustersReads[it->first];
            } else{
                //这种的情况是ISM类, 不用管是哪种ISM; //如果后续改识别ISM类型的时候, 更改这里; ★★★★★★★
                FSMISMO.ISM[it->first] = ClustersReads[it->first];
                //这些ISM是由哪些FSM得来的, 有的FSM可能已经有number, 有的没有number的就记为0;
                for (int i = 0; i < statis_ISM.size(); i++){
                    auto isIn = FSMISMO.FSM.find(statis_ISM[i]);
                    if (isIn == FSMISMO.FSM.end()){
                        FSMISMO.FSM[statis_ISM[i]].first = GroupAnno[statis_ISM[i]];
                        FSMISMO.FSM[statis_ISM[i]].second = {};
                    }
                }
            }
        }
    }
    return FSMISMO;
}



struct HighLowClusters{
    //所有的需要进一步进行识别的Cluster, 以及他们对应的数量;
    std::map<size_t, std::vector<std::string>> HighConClusters;
    //淘汰的Cluster及其他们对应的数量;
    std::map<size_t, std::vector<std::string>> LowConClusters;
    std::unordered_map<size_t, std::string> HighStrand;
};

HighLowClusters get_HighLow_clusters(std::unordered_map<std::size_t, std::vector<std::string>>& Others_Clusters, 
                                    std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs, 
                                    std::unordered_map<std::string, std::vector<int>>& readsigns,
                                    const int& sjsupportnumbers){
    //结构体初始化;
    HighLowClusters HighLowC = {};
    std::map<std::size_t, std::vector<std::string>> onlyHighC;
    std::unordered_map<size_t, std::string> onlyHighStrand;
    //在for循环的时候要删除;
    std::map<std::size_t, std::vector<std::string>> onlyHighC_copy;
    std::unordered_map<size_t, std::string> onlyHighStrand_copy;

    //读取的变量;
    std::vector<std::string> each_cluster; //每个cluster的read的名称;
    std::size_t hashname; //每个cluster的名称;
    std::vector<int> eachreadsign;
    
    //中间用到的变量;
    int SumCount = 0;
    int StrandCount = 0;
    int pos1 = 0;
    int pos2 = 0;
    std::vector<int> results;
    int SjNumber = 0;
    bool StrandFlag = 0;

    for (auto it = Others_Clusters.begin(); it != Others_Clusters.end(); ++it) {
        hashname = it->first; //Cluster对应的哈希值;
        each_cluster = it->second; //Cluster对应的所有read名称的vector;
        //首先判断cluster中reads的数量, 一个的很多;
        if (each_cluster.size() > 1){
            // 多于1个的继续考虑;
            SumCount = 0;
            StrandCount = 0;
            SjNumber = readsigns[each_cluster[0]].size();
            results.assign(SjNumber, 0);
            StrandFlag = 0;

            for (const auto& each:each_cluster){
                eachreadsign = readsigns[each];
                for (int i = 0; i < eachreadsign.size(); i++) {
                    if (eachreadsign[i] == 1) {
                        results[i] = results[i] + eachreadsign[i];
                        if (StrandCount == 2) {
                            StrandFlag = 1;
                        }
                        StrandCount = 1;
                    } else if (eachreadsign[i] == 2) {
                        results[i] = results[i] + 1;
                        if (StrandCount == 1) {
                            StrandFlag = 1;
                        }
                        StrandCount = 2;
                    }
                }
            }

            for (const auto& eachBit:results) {
                if (eachBit > sjsupportnumbers) {
                    SumCount = SumCount + 1;
                }
            }

            // 看看满足条件;
            if (SumCount == SjNumber) {
                if (StrandFlag) {
                    // 链不一致;
                    HighLowC.LowConClusters[hashname] = Others_Clusters[hashname];               
                } else {
                    if (StrandCount == 1) {
                        onlyHighStrand[hashname] = "+";
                        onlyHighC[hashname] = Others_Clusters[hashname];
                    } else {
                        onlyHighStrand[hashname] = "-";
                        onlyHighC[hashname] = Others_Clusters[hashname];
                    }
                }
            } else {
                HighLowC.LowConClusters[hashname] = Others_Clusters[hashname];
                // std::cout << "Low类不满足>2" << std::endl;
            }
        } else {
            //一个类里面如果等于1个直接就不需要考虑了, 保存到Low的里面, 后面等回收;
            HighLowC.LowConClusters[hashname] = Others_Clusters[hashname];
            // std::cout << "Low类小于1" << std::endl;
        }
    }
    // std::cout << onlyHighC.size() << " High个  ***  " << HighLowC.LowConClusters.size() << " Low个" << std::endl;
    onlyHighC_copy = onlyHighC;
    onlyHighStrand_copy = onlyHighStrand;
    //还有一步, 将短的cluster去掉, 留下长的cluster;
    std::vector<std::array<int,2>> SJs1;
    std::vector<std::array<int,2>> SJs2;
    for (const auto& clus1: onlyHighC) {
        SJs1 = AllSJs[clus1.second[0]];
        for (const auto& clus2: onlyHighC){
            SJs2 = AllSJs[clus2.second[0]];
            if (SJs2.size()>SJs1.size()) {
                if (std::includes(SJs2.begin(), SJs2.end(), SJs1.begin(), SJs1.end())){                  
                    //这种情况是clus1是clus2的子集;
                    //判断成功连续就可以了;
                    pos1 = std::find(SJs2.begin(), SJs2.end(), SJs1[0]) - SJs2.begin();
                    pos2 = std::find(SJs2.begin(), SJs2.end(), SJs1[SJs1.size()-1]) - SJs2.begin();
                    if (pos2-pos1 == SJs1.size()-1){
                        //这种情况说明1是2的子集;
                        HighLowC.LowConClusters[clus1.first] = clus1.second;
                        onlyHighC_copy.erase(clus1.first);
                        onlyHighStrand_copy.erase(clus1.first);
                    }
                }
            }
        }
    }

    for (const auto& eachHighCluster:onlyHighC_copy){
        HighLowC.HighConClusters[eachHighCluster.first] = eachHighCluster.second;
        HighLowC.HighStrand[eachHighCluster.first] = onlyHighStrand_copy[eachHighCluster.first];
    }
    return HighLowC;
}



int IntervalIntersection(std::vector<std::array<int,2>> firstList, std::vector<std::array<int,2>> secondList) {
    std::vector<std::array<int,2>> res; //定义返回的结果;
    int result = 0;
    int i = 0;
    int j = 0;

    while(i < firstList.size() && j < secondList.size()){
        int start = std::max(firstList[i][0], secondList[j][0]);
        int end = std::min(firstList[i][1], secondList[j][1]);
        if(start <= end){
            res.push_back({start, end});
            result = result + (end-start+1);       
        }
        if(firstList[i][1] < secondList[j][1]){
            i++;
        }else{
            j++;
        }
    }
    return result;
}


int IntervalMerge(std::vector<std::array<int,2>> firstList, std::vector<std::array<int,2>> secondList) {
    //将两个序列集合合并为一个;
    std::vector<std::array<int,2>> intervals;
    intervals.resize(firstList.size()+secondList.size());
    merge(firstList.begin(), firstList.end(), secondList.begin(), secondList.end(), intervals.begin());
    
    sort(intervals.begin(), intervals.end()); //对所有的区间从小到大排序;
    std::vector<std::array<int,2>> merged; //初始化输出结果;
    int result = 0;
    
    for (int i=0; i<intervals.size();) {  
        int left = intervals[i][0];
        int right = intervals[i][1];
        while (i+1 < intervals.size() && right >= intervals[i+1][0]){
            right = std::max(right, intervals[i+1][1]);
            i++;
        }
        merged.push_back({left, right});
        result = result + (right-left+1);
        i++;
    }
    // std::cout << "结果是:" << result << std::endl;
    return result;
}


void get_filtered_FSM(std::map<std::size_t, std::vector<std::string>>& LowReads, 
                      std::unordered_map<std::string, std::vector<std::array<int,2>>>& GroupAnno, 
                      std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>>& AllFSM, 
                      std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs, 
                      std::unordered_map<std::string, std::string>& GroupReadFile,
                      std::ofstream& Trace){
    // std::mutex& TraceMutex
    //需要保存的变量;
    std::unordered_map<std::string, std::unordered_map<std::size_t, std::vector<std::string>>> RecycleAnnoName_value;
    std::set<std::size_t> DeleteItem;
    //每次计算出相差的距离;
    int RecycleLength = 0;
    int count = 0;
    std::string minAnnoName;
    int dis = 10000;
    //临时变量, 保存每个cluster潜在的多个FSM;
    std::unordered_map<std::string, int> potential_FSM;
    std::vector<std::array<int,2>> LowReadSJVec;

    std::string first_part;
    std::string second_part;

    //需要循环每一个过滤的cluster, 每个cluster可能有几个潜在的FSM;
    for (const auto& each_Readcluster:LowReads){
        LowReadSJVec = AllSJs[each_Readcluster.second[0]];
        // std::cout << each_anno.first << std::endl;
        //循环所有的过滤出去的reads_cluster;
        potential_FSM.clear();
        for (const auto& each_anno:GroupAnno){
            //如果注释与reads的cluster相同长度才可以计算, 为了减少计算的数量;
            if (each_anno.second.size() == LowReadSJVec.size()){
                RecycleLength = IntervalMerge(each_anno.second, LowReadSJVec) - IntervalIntersection(each_anno.second, LowReadSJVec);
                if (RecycleLength < 15*LowReadSJVec.size()){
                    //如果距离符合要求;
                    potential_FSM[each_anno.first] = RecycleLength;
                }
            }
        }
        //对于每个cluster查看有几个潜在的FSM;
        //如果只有1个是最好办的情况;
        if (potential_FSM.size() == 1){
            //首先将这个cluster标记;
            DeleteItem.insert(each_Readcluster.first);
            //输出到追踪文件;
            size_t pos = (*(potential_FSM.begin())).first.find('|');
            if (pos != std::string::npos){
                first_part = (*(potential_FSM.begin())).first.substr(0, pos);
                second_part = (*(potential_FSM.begin())).first.substr(pos + 1);
            }
            if (Trace.is_open()) {
                for (const auto& EachRead:each_Readcluster.second){
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << EachRead << '\t' << "FSM" << '\t' << second_part << '\t' << first_part << '\t' << GroupReadFile[EachRead] << '\n'; 
                }
            }                        

            //查看这个潜在的FSM是不是已经出现过;
            if (RecycleAnnoName_value.find((*potential_FSM.begin()).first) != RecycleAnnoName_value.end()){
                //如果已经出现过, 只需要在原来的基础上加上数值;
                RecycleAnnoName_value[(*potential_FSM.begin()).first][each_Readcluster.first] = each_Readcluster.second;
            } else{
                //如果没有出现过;
                RecycleAnnoName_value[(*potential_FSM.begin()).first] = {};
                RecycleAnnoName_value[(*potential_FSM.begin()).first][each_Readcluster.first] = each_Readcluster.second;
            }
        //如果有多个, 则取距离最小的
        } else if (potential_FSM.size() > 1) {
            //首先将这个cluster标记;
            DeleteItem.insert(each_Readcluster.first);
            //找到距离最小的;
            dis = 100000;
            for (const auto& Po:potential_FSM){
                if (Po.second < dis){
                    dis = Po.second;
                    minAnnoName = Po.first;
                }
            }
            //输出追踪文件;
            size_t pos = minAnnoName.find('|');
            if (pos != std::string::npos){
                first_part = minAnnoName.substr(0, pos);
                second_part = minAnnoName.substr(pos + 1);
            }            
            if (Trace.is_open()) {
                for (const auto& EachRead:each_Readcluster.second){
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << EachRead << '\t' << "FSM" << '\t' << second_part << '\t' << first_part << '\t' << GroupReadFile[EachRead] << '\n'; 
                }
            }   

            //查看这个潜在的FSM——minAnnoName是不是已经出现过?
            if (RecycleAnnoName_value.find(minAnnoName) != RecycleAnnoName_value.end()){
                //如果已经出现过, 只需在原来的基础上加上数值;
                RecycleAnnoName_value[minAnnoName][each_Readcluster.first] = each_Readcluster.second;
            } else{
                //如果没有出现过;
                RecycleAnnoName_value[minAnnoName] = {};
                RecycleAnnoName_value[minAnnoName][each_Readcluster.first] = each_Readcluster.second;
            }
        }
        //整个cluster完毕;
    }   
    // int new_count = 0;
    std::vector<std::string> newName;

    for (const auto& newFSM:RecycleAnnoName_value){
        //如果之前的FSM中已经存在;
        // new_count = 0;
        newName.clear();
        for (const auto& allCluster:newFSM.second){ //每个类;
            // new_count = new_count + allCluster.second.size();
            newName.insert(newName.end(), allCluster.second.begin(), allCluster.second.end());
        }

        if (AllFSM.find(newFSM.first) != AllFSM.end()){

            AllFSM[newFSM.first].second.insert(AllFSM[newFSM.first].second.end(), newName.begin(), newName.end());
        } else{

            AllFSM[newFSM.first].first = GroupAnno[newFSM.first];

            AllFSM[newFSM.first].second.insert(AllFSM[newFSM.first].second.end(), newName.begin(), newName.end());
        }
    }

    for (const auto& eachSJ:DeleteItem){
        auto it = LowReads.find(eachSJ);
        if (it != LowReads.end()){
            LowReads.erase(it);
        } else { //后期可以去掉;
            std::cout << "*** That doesn't make sense. There should be a Bug! ***" << std::endl;
        }
    }
}




struct SpliceChainClass
{
    std::unordered_map<std::size_t, std::array<int,2>> ClusterCoverage;
    std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>> FSM;
    std::map<std::size_t, std::vector<std::string>> ISM;
    std::map<size_t, std::vector<std::string>> HighConClusters;
    std::unordered_map<size_t, std::string> HighStrand;
};

SpliceChainClass generate_splice_chain_class(const int& dataMode, 
                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs, 
                                std::unordered_map<std::string, std::array<int,2>>& groupreadcoverage, 
                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupannotations, 
                                std::unordered_map<std::string, std::vector<int>> groupreadsigns,
                                std::unordered_map<std::string, std::array<int,2>>& AnnoCoverage,
                                std::unordered_map<std::string, std::string>& groupreadfiles,
                                std::ofstream& traceFilePath,
                                const int& Sj_Support_Number,
                                std::string& chrname) {
    // std::mutex& readTraceMutex;
    SpliceChainClass SCC;
    ReferenceCluster FsmIsmOthers;
    HighLowClusters Others2HighLow;
    std::unordered_map<std::size_t, std::vector<std::string>> groupCluster = classifyReadsVec(groupreadsjs);
    // FSM ISM 根本没必要校正区间; 真正需要校正的是Low的;
    // 后面有一步回收, 能不能达到相同效果? 
    SCC.ClusterCoverage = get_every_cluster_begin_end(groupCluster, groupreadcoverage);
    FsmIsmOthers = get_FSM_and_others(groupCluster, groupannotations, AnnoCoverage, groupreadcoverage, groupreadsjs, groupreadfiles, traceFilePath); //readTraceMutex;
    Others2HighLow = get_HighLow_clusters(FsmIsmOthers.Others, groupreadsjs, groupreadsigns, Sj_Support_Number);
     
    get_filtered_FSM(Others2HighLow.LowConClusters, groupannotations, FsmIsmOthers.FSM, groupreadsjs, groupreadfiles, traceFilePath); //readTraceMutex;     
    
    // 将有用结果输出;
    SCC.FSM = FsmIsmOthers.FSM;
    SCC.ISM = FsmIsmOthers.ISM;
    SCC.HighConClusters = Others2HighLow.HighConClusters;
    SCC.HighStrand = Others2HighLow.HighStrand;

    return SCC;
}


int findNonOverlappingIntervals(const std::vector<std::array<int, 2>>& intervals1, const std::vector<std::array<int, 2>>& intervals2) {
    std::vector<std::array<int, 2>> result;
    int finally_number = 0;

    for (const auto& interval1 : intervals1) {
        bool hasOverlap = false;
        for (const auto& interval2 : intervals2) {
            if (!(interval1[1] < interval2[0] || interval1[0] > interval2[1])) {
                // Overlapping, don't include this interval
                hasOverlap = true;
                break;
            }
        }
        if (!hasOverlap) {
            result.push_back(interval1);
            finally_number = finally_number + (interval1[1]-interval1[0]+1);
        }
    }
    for (const auto& interval2 : intervals2) {
        bool hasOverlap = false;

        for (const auto& interval1 : intervals1) {
            if (!(interval2[1] < interval1[0] || interval2[0] > interval1[1])) {
                // Overlapping, don't include this interval
                hasOverlap = true;
                break;
            }
        }
        if (!hasOverlap) {
            result.push_back(interval2);
            finally_number = finally_number + (interval2[1]-interval2[0]+1);
        }
    }
    return finally_number;
}

struct DistanceInform
{
    std::unordered_map<int, std::set<int>> NodeDistanceSet;
    //序号名称和节点名称的字典;
    std::unordered_map<int, std::vector<std::array<int,2>>> Index2Unknown;
    std::unordered_map<int, std::string> Index2Anno;
    std::unordered_map<int, int> Index2Count;
    std::unordered_map<int, std::size_t> Index2hashname;
    std::unordered_map<int, std::string> Index2novelname; //后面才用到;
};


DistanceInform get_distance_matrix(std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupannotations, 
                                   std::map<std::size_t, std::vector<std::string>>& HighClusters, 
                                   std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs, int DistanceG, 
                                   std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>>& AllFSM){
    //输出结果的二维矩阵;
    DistanceInform Graph_Init_inform;
    //初始化
    for (int i = 0; i < groupannotations.size()+HighClusters.size(); ++i){
        Graph_Init_inform.NodeDistanceSet[i] = {};
    }

    std::vector<std::array<int,2>> SJVec1;
    std::vector<std::array<int,2>> SJVec2;
    int LineElement = 0; //距离;
    int NameNumber = 0;
    int HighNumber = groupannotations.size();

    //首先计算的一部分是注释的isoform与候选的isoform之间的距离; 相当于整个矩阵的part2;
    //如果注释这部分没有;
    if (groupannotations.size() != 0){
        for (const auto& EachAnno : groupannotations){
            Graph_Init_inform.Index2Anno[NameNumber] = EachAnno.first;
            //提取这个注释是否在FSM中;
            auto it = AllFSM.find(EachAnno.first);
            //检查是否找到,如果在FSM中, 则将数量提取, 否则, 数量记为0;
            if (it != AllFSM.end()) {
                Graph_Init_inform.Index2Count[NameNumber] = (it->second).second.size();
            } else{
                Graph_Init_inform.Index2Count[NameNumber] = 0;
            }
            //清空每一行的距离重新赋值;
            if (HighClusters.size() != 0){
                //初始化;
                HighNumber = groupannotations.size();
                for (const auto& eachHC : HighClusters){
                    //从这里开始计算两两之间的距离大小;
                    SJVec1 = groupreadsjs[eachHC.second[0]];
                    LineElement = IntervalMerge(EachAnno.second, SJVec1) - IntervalIntersection(EachAnno.second, SJVec1) - findNonOverlappingIntervals(EachAnno.second, SJVec1);
                    // std::cout << LineElement << std::endl;
                    if ((LineElement < DistanceG+1) && (LineElement > 0)){
                        Graph_Init_inform.NodeDistanceSet[NameNumber].insert(HighNumber);
                        Graph_Init_inform.NodeDistanceSet[HighNumber].insert(NameNumber);
                    }
                    HighNumber++;
                }
            }
            NameNumber++;
        }
    }

    if (HighClusters.size() != 0) {
        //接下来计算的候选的isoform之间的距离;通过part2得到part24;
        for (const auto& HC1 : HighClusters){
            SJVec1 = groupreadsjs[HC1.second[0]];
            Graph_Init_inform.Index2Unknown[NameNumber] = SJVec1;
            Graph_Init_inform.Index2Count[NameNumber] = HC1.second.size();
            Graph_Init_inform.Index2hashname[NameNumber] = HC1.first;
            HighNumber = groupannotations.size();
            for (const auto& HC2 : HighClusters){
                SJVec2 = groupreadsjs[HC2.second[0]];
                LineElement = IntervalMerge(SJVec1, SJVec2) - IntervalIntersection(SJVec1, SJVec2) - findNonOverlappingIntervals(SJVec1, SJVec2);
                if ((LineElement < DistanceG+1) && (LineElement > 0)){
                    Graph_Init_inform.NodeDistanceSet[NameNumber].insert(HighNumber);
                }
                HighNumber++;
            }
            NameNumber++;
        }
    }
    return Graph_Init_inform;
}


int find_max_neighborhood_node(std::set<int>& fireNode, std::unordered_map<int, std::set<int>>& evidence){
    int node_max_nei = 0;
    int max_number = 0;
    for (const auto& node:fireNode){
        if (evidence[node].size() > max_number){
            node_max_nei = node;
            max_number = evidence[node].size();
        }
    }
    return node_max_nei;
}


std::vector<std::vector<int>> Find_Maximal_Cliques(std::unordered_map<int, std::set<int>>& adj){
    std::vector<std::vector<int>> CliqueVec;
    //初始化Q;
    std::vector<int> Q = {};
    std::set<int> cand;
    std::set<int> subg;
    std::set<int> ext_u;
    std::vector<std::tuple<std::set<int>, std::set<int>, std::set<int>>> stack;

    int q = 0;
    int max_node = 0;
    std::set<int> adj_q;
    std::set<int> subg_q;
    std::set<int> cand_q;

    if (adj.size() == 0){
        return CliqueVec;
    } else{
        //adj已经有了;
        for (const auto& eachnode:adj) cand.insert(eachnode.first);
        subg = cand;
        Q.push_back(-1);
        //首先找到邻居节点最多的节点;        
        max_node = find_max_neighborhood_node(subg, adj);
        //使用std::set_difference算法计算差集;
        std::set_difference(cand.begin(), cand.end(),
                            adj[max_node].begin(), adj[max_node].end(),
                            std::inserter(ext_u, ext_u.begin()));

        while (true)
        {
            if (!ext_u.empty()) 
            {
                q = *ext_u.begin(); //获取集合中的第一个元素
                ext_u.erase(q);
                cand.erase(q);
                Q.back() = q;
                adj_q = adj[q];
                subg_q.clear();
                std::set_intersection(subg.begin(), subg.end(),
                                    adj_q.begin(), adj_q.end(),
                                    std::inserter(subg_q, subg_q.begin())); //subg_q
                if (subg_q.empty())
                {
                    CliqueVec.push_back(Q);
                } 
                else
                {   
                    cand_q.clear();
                    std::set_intersection(cand.begin(), cand.end(),
                                    adj_q.begin(), adj_q.end(),
                                    std::inserter(cand_q, cand_q.begin())); //cand_q
                
                    if (!cand_q.empty())
                    {
                        stack.push_back(std::make_tuple(subg, cand, ext_u));
                        Q.push_back(-1); //添加一个占位符;
                        subg = subg_q;
                        cand = cand_q;
                        max_node = find_max_neighborhood_node(subg, adj);
                        ext_u.clear();
                        std::set_difference(cand.begin(), cand.end(),
                            adj[max_node].begin(), adj[max_node].end(),
                            std::inserter(ext_u, ext_u.begin()));
                    } 
                }
            }
            else
            {   
                if (stack.size() != 0 && Q.size() != 0)
                {
                    Q.pop_back();
                    std::tie(subg, cand, ext_u) = stack.back();
                    stack.pop_back();
                } 
                else
                {
                    break;
                }
            }
        }
        return CliqueVec;
    }
}



bool compareBySize(const std::vector<int>& a, const std::vector<int>& b) {
    // std::cout << "Comparing sizes: " << a.size() << " and " << b.size() << std::endl;
    if (a.empty() && b.empty()) {
        return false; // 两个向量都为空，返回相等或任何你需要的值
    }
    if (a.empty()) {
        return true; // 空向量被认为比非空向量小
    }
    if (b.empty()) {
        return false; // 非空向量被认为比空向量大
    }
    auto AminElement = std::min_element(a.begin(), a.end());
    auto BminElement = std::min_element(b.begin(), b.end());
    // 比较最小元素
    if (*AminElement != *BminElement) {
        return *AminElement < *BminElement; // 返回最小元素的比较结果
    }
    // 如果最小元素相等，按大小比较
    return a.size() < b.size();
}



struct DetectionResults
{
    std::set<int> NodesKnownTrue;
    std::set<int> TrueNodeSet;
    std::set<int> FalseNodeSet;
};


DetectionResults Transcript_Detection(std::vector<std::vector<int>>& CqVec, 
                                      std::unordered_map<int, std::string>& Number2Anno, 
                                      std::unordered_map<int, std::vector<std::array<int,2>>>& Number2Unknown, 
                                      std::unordered_map<int, int>& Number2Count) {
    //整个染色体中已知正确的节点;
    //定义结果, 识别为正确的和识别为错误的节点;
    //初始化结构体;
    DetectionResults Node_Results = {};

    //函数用到的中间变量;
    std::set<int> subgraph_known_true_nodes;
    std::set<int> subgraph_remain_nodes;
    int standard_read_number = std::numeric_limits<int>::max();
    int known_true_number = 0;
    int max_node_value = 0;


    //遍历一遍已知注释的点, 将其存入NodesKnownTrue的Set中;
    for (const auto& Anno : Number2Anno){
        Node_Results.NodesKnownTrue.insert(Anno.first);
    }

    if (CqVec.size() != 0){

        if (CqVec.size() > 1) {
            std::sort(CqVec.begin(), CqVec.end(), compareBySize);
        }
        //遍历每个团;
        for (const auto& eachClique : CqVec){
            standard_read_number = std::numeric_limits<int>::max();

            std::set<int> CliqueSet(eachClique.begin(), eachClique.end());
            subgraph_known_true_nodes.clear();
            std::set_intersection(CliqueSet.begin(), CliqueSet.end(), Node_Results.NodesKnownTrue.begin(), Node_Results.NodesKnownTrue.end(),
                          std::inserter(subgraph_known_true_nodes, subgraph_known_true_nodes.begin()));
            subgraph_remain_nodes.clear();
            std::set_difference(CliqueSet.begin(), CliqueSet.end(), subgraph_known_true_nodes.begin(), subgraph_known_true_nodes.end(),
                        std::inserter(subgraph_remain_nodes, subgraph_remain_nodes.begin()));

            if (subgraph_known_true_nodes.size() != 0){

                if (subgraph_known_true_nodes.size() == 1){
                    standard_read_number = Number2Count[(*subgraph_known_true_nodes.begin())];
                } else {
                    //如果有多个正确的节点;
                    for(const auto& each_True_node : subgraph_known_true_nodes){
                        if (Number2Count[each_True_node] < standard_read_number) {
                            standard_read_number = Number2Count[each_True_node];
                        }
                    }
                }

                if (subgraph_remain_nodes.size() != 0){ 
                    //有已知相邻的节点, 即使number为0也可以加入识别的;
                    for (const auto& eachknownnode:subgraph_known_true_nodes){
                        Node_Results.TrueNodeSet.insert(subgraph_known_true_nodes.begin(), subgraph_known_true_nodes.end());
                    }
                    //利用已知节点, 看一下其他的节点, 不然直接剔除;
                    for (const auto& each_remain_node : subgraph_remain_nodes){
                        //查看是否满足条件;
                        if (standard_read_number == 0){
                            Node_Results.FalseNodeSet.insert(each_remain_node);
                        } else {
                            if (Number2Count[each_remain_node] < standard_read_number){
                                Node_Results.FalseNodeSet.insert(each_remain_node);
                            } else {
                                Node_Results.TrueNodeSet.insert(each_remain_node);
                            }
                        }
                    }
                } else { 

                    for (const auto& eachknownnode:subgraph_known_true_nodes){
                        if (Number2Count[eachknownnode] > 0){
                            Node_Results.TrueNodeSet.insert(eachknownnode);
                        }
                    } 
                }

            } else {

                if (subgraph_remain_nodes.size() == 1){
                    //只有一个节点, 暂时认为对;
                    if (Number2Count[*subgraph_remain_nodes.begin()] > 0) {
                        Node_Results.TrueNodeSet.insert(*subgraph_remain_nodes.begin());
                    }
                } else{

                    max_node_value = 0;
                    for (const auto& each_remain_node : subgraph_remain_nodes){
                        if (max_node_value < Number2Count[each_remain_node]){
                            max_node_value = Number2Count[each_remain_node];
                        }
                    }
                    for (const auto& each_remain_node : subgraph_remain_nodes){
                        if (max_node_value == Number2Count[each_remain_node]){
                            Node_Results.TrueNodeSet.insert(each_remain_node);
                        }
                        else
                        {
                            Node_Results.FalseNodeSet.insert(each_remain_node);
                        }
                    }   
                }
            }
            //每个子图一次循环;
        }

        for (const auto& eachFalse:Node_Results.FalseNodeSet){
            auto iter = Node_Results.TrueNodeSet.find(eachFalse);
            if (iter != Node_Results.TrueNodeSet.end()){
                Node_Results.TrueNodeSet.erase(eachFalse);
            }
        }
    }
    return Node_Results;
}


struct Solvent {
    std::map<std::string, std::unordered_map<std::string, int>> File_FSM;
    std::map<std::string, std::map<std::size_t, std::vector<std::string>>> File_ISM;
    std::map<std::string, std::map<std::size_t, std::vector<std::string>>> File_HighConClusters;
};

Solvent get_Solvent_FsmIsmHigh(SpliceChainClass& FsmIsmHigh, int& FileNumber, GroupInformation& groupinformations) {
    Solvent FileSpliceChains;
    std::unordered_map<std::string, int> thisTranscript_File2Count; // FSM的临时变量;
    std::map<std::string, std::vector<std::string>> thisTranscript_File2Count_ISM; // ISM和High的临时变量;

    // 拆分FSM;
    if (FsmIsmHigh.FSM.size() > 0) {
        // 首先初始化;
        for (int i = 0; i < FileNumber; i++) {
            FileSpliceChains.File_FSM[std::to_string(i)] = {};
        }        
        // 对每个FSM循环;
        for (const auto& eachCluster:FsmIsmHigh.FSM) {
            // std::cout << "总 " << eachCluster.first << ": " << eachCluster.second.second.size() << std::endl;
            thisTranscript_File2Count.clear();
            for (int i = 0; i < FileNumber; i++) {
                thisTranscript_File2Count[std::to_string(i)] = 0;
            }
            std::string AnnoName = eachCluster.first;
            // 对这个transcript所属的每条reads进行循环;
            for (const auto& eachRead:eachCluster.second.second) {
                // 这条read所在的文件;
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                thisTranscript_File2Count[read_file] = thisTranscript_File2Count[read_file] + 1;
            }
            // 所有的reads循环结束;
            // 结果填写到FSM中;
            for (const auto& eachFile:thisTranscript_File2Count) {
                // std::cout << "File: " << eachFile.first << " " << eachFile.second << std::endl;
                if (eachFile.second > 0) {
                    FileSpliceChains.File_FSM[eachFile.first][AnnoName] = eachFile.second;
                }
            }            
        }
    } // 结束;

    // 拆分ISM;
    if (FsmIsmHigh.ISM.size() > 0) {
        // 首先初始化;
        for (int i = 0; i < FileNumber; i++) {
            FileSpliceChains.File_ISM[std::to_string(i)] = {};
        }
        // 对每个ISM循环;
        for (const auto& eachCluster:FsmIsmHigh.ISM) {
            thisTranscript_File2Count_ISM.clear();
            for (int i = 0; i < FileNumber; i++) {
                thisTranscript_File2Count_ISM[std::to_string(i)] = {};
            }
            std::size_t NonName = eachCluster.first;
            // 对这个ISM所属的每条reads进行循环;
            for (const auto& eachRead:eachCluster.second) {
                // 这条read所在的文件;
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                thisTranscript_File2Count_ISM[read_file].push_back(eachRead);
            }
            // 所有的reads循环结束;
            // 结果填写到ISM中;
            for (const auto& eachFile:thisTranscript_File2Count_ISM) {
                if (eachFile.second.size() > 0) {
                    FileSpliceChains.File_ISM[eachFile.first][NonName] = eachFile.second;
                }
            }
        }
    }

    // 拆分HighConfidenceClusters;
    if (FsmIsmHigh.HighConClusters.size() > 0) {
        // 首先初始化;
        for (int i = 0; i < FileNumber; i++) {
            FileSpliceChains.File_HighConClusters[std::to_string(i)] = {};
        }
        // 对每个HighConCluster循环;
        for (const auto& eachCluster:FsmIsmHigh.HighConClusters) {
            thisTranscript_File2Count_ISM.clear();
            for (int i = 0; i < FileNumber; i++) {
                thisTranscript_File2Count_ISM[std::to_string(i)] = {};
            }
            std::size_t NonName = eachCluster.first;
            // 对这个HighConCluster所属的每条reads进行循环;
            for (const auto& eachRead:eachCluster.second) {
                // 这条read所在的文件;
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                thisTranscript_File2Count_ISM[read_file].push_back(eachRead);
            }
            // 所有的reads循环结束;
            // 结果填写到HighConCluster中;
            for (const auto& eachFile:thisTranscript_File2Count_ISM) {
                if (eachFile.second.size() > 0) {
                    FileSpliceChains.File_HighConClusters[eachFile.first][NonName] = eachFile.second;
                }
            }             
        }
    }
    return FileSpliceChains;
}


std::string get_co_sj_gene(std::vector<std::array<int,2>>& novel_isoform, 
                            std::set<std::string>& gene_set,
                            std::map<std::string, std::vector<std::string>>& gene2tx,
                            std::unordered_map<std::string, std::vector<std::array<int,2>>>& tx2sj) {
    //
    std::vector<std::string> TempGene;
    std::vector<int> TempDist;
    std::set<std::array<int,2>> SJ_set; 
    std::vector<std::string> ALL_tx;
    std::set<std::array<int, 2>> novel_isoform_set(novel_isoform.begin(), novel_isoform.end());
    std::set<std::array<int, 2>> intersection;

    // 对每个基因循环;
    for (const auto& A_Gene:gene_set) {
        SJ_set.clear();
        intersection.clear();
        TempGene.push_back(A_Gene);
        ALL_tx = gene2tx[A_Gene];
        for (const auto& A_Tx:ALL_tx) {
            // 这个转录本在里面;
            if (tx2sj.find(A_Tx) != tx2sj.end()) {
                SJ_set.insert(tx2sj[A_Tx].begin(), tx2sj[A_Tx].end());
            }
        }
        std::set_intersection(novel_isoform_set.begin(), novel_isoform_set.end(), SJ_set.begin(), SJ_set.end(),
                          std::inserter(intersection, intersection.begin()));
        TempDist.push_back(intersection.size());
    }
    // 判断哪个基因多;
    auto max_it = std::max_element(TempDist.begin(), TempDist.end());
    std::string MaxGeneName;
    if (max_it != TempDist.end()) {
        if (*max_it > 0) {
            MaxGeneName = TempGene[std::distance(TempDist.begin(), max_it)];
        } else {
            MaxGeneName = "";
        }
    } else {
        MaxGeneName = "";
    }    
    return MaxGeneName;
}


struct OutputInformation
{
    std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, double>> Transcript_Annotations;
    std::unordered_map<std::string, std::string> transcript2gene;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> File_TranscriptNumber;
};

//18. 将识别的结果输出到txt;
OutputInformation Write_Detection_Transcript2gtf(std::ofstream& Updated_Files, std::ofstream& Trace,
                                                DetectionResults& noderesults, DistanceInform& Disinform, 
                                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupannotations, 
                                                std::map<std::string, std::vector<std::array<int,2>>>& Annoexon, 
                                                std::unordered_map<std::size_t, std::array<int, 2>>& groupreadcoverage, 
                                                std::set<std::string>& groupallgene, 
                                                std::map<std::string, std::array<int,2>>& Annogenecovergae, 
                                                std::unordered_map<std::string, std::string>& GTF_Transcript_Strand, 
                                                std::unordered_map<size_t, std::string>& High_Strand, 
                                                std::string& chrname, std::string& group_size, 
                                                std::map<size_t, std::vector<std::string>>& HighClusters,
                                                std::unordered_map<std::string, std::string>& groupreadfiles,
                                                std::unordered_map<std::string, std::array<int,2>>& SingleExonAnno,
                                                std::map<std::string, std::vector<std::string>>& gtf_gene2tx) {
    // std::mutex& GtfMutex, std::mutex& TraceMutex;
    OutputInformation FinalAnnotations;
    std::string itsname;
    std::vector<std::array<int,2>> itssj;
    std::vector<std::array<int,2>> itsexon;
    std::vector<std::string> TempGene;
    std::string TempGeneName;
    std::vector<int> TempDist;
    int novel_count = 0;
    std::size_t nameya;
    std::array<int,2> BE;
    std::string first_part;
    std::string second_part;

    if (Updated_Files.is_open()){
        // 第一次
        for (const auto aaa:noderesults.TrueNodeSet){
            if (aaa < Disinform.Index2Anno.size()){
                //有注释的部分;
                itsname = Disinform.Index2Anno[aaa];
                itssj = groupannotations[itsname];
                itsexon = Annoexon[itsname];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];
                
                size_t pos = itsname.find('|');
                if (pos != std::string::npos) {
                    // 根据 "|" 将字符串分割为两部分
                    first_part = itsname.substr(0, pos);
                    second_part = itsname.substr(pos + 1);
                }
                {
                    std::unique_lock<std::mutex> lock(updatedGtfMutex);
                    Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t' << itsexon[0][0] << '\t' << itsexon[itsexon.size()-1][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\";" << '\n';
                    for (int i = 0; i < itsexon.size(); i++){
                        Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t' << itsexon[i][0] << '\t' << itsexon[i][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\"; exon_number \"" << i << "\";" << '\n';
                    }
                }
                FinalAnnotations.transcript2gene[itsname] = first_part;

            } else {
                novel_count = novel_count + 1;
                //没有注释的部分;
                itsname = chrname + "-novel-" + group_size + "-" + std::to_string(aaa) + "-" + std::to_string(novel_count);
                Disinform.Index2novelname[aaa] = itsname;
                itssj = Disinform.Index2Unknown[aaa];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];
                nameya = Disinform.Index2hashname[aaa];
                BE = groupreadcoverage[nameya];

                TempGene.clear();
                // gene set里面已经有很多isoform; 当前isoform的sj就是itssj;
                if (groupallgene.size() != 0) {
                    TempGeneName = get_co_sj_gene(itssj, groupallgene, gtf_gene2tx, groupannotations);

                    if (TempGeneName.size() != 0) {
                        first_part = TempGeneName;
                    } else {
                        // 这个时候就是多个exon的注释和单个exon的注释, 从里面选一个;
                        for (const auto& every_gene:groupallgene) {
                            if (Annogenecovergae[every_gene][0] <= BE[1] && BE[0] <= Annogenecovergae[every_gene][1]) {
                                TempGene.push_back(every_gene);
                            }
                        }
                        //novel isoform没有候选的基因;
                        if (TempGene.size() == 0) {
                            first_part = "NA";
                            // novel isoform只有1个候选的基因;
                        } else if (TempGene.size() == 1) {
                            first_part = TempGene[0];
                            //novel isoform有2个以上候选的基因;
                        } else {
                            TempDist.clear();
                            for (const auto& each_gene:TempGene){
                                int a = abs(BE[1] - Annogenecovergae[each_gene][1]);
                                int b = abs(Annogenecovergae[each_gene][0] - BE[0]);
                                TempDist.push_back((a>b)?a:b);
                            }   
                            auto min_it = std::min_element(TempDist.begin(), TempDist.end());
                            first_part = TempGene[std::distance(TempDist.begin(), min_it)];                           
                        }
                    }

                } else {
                    // 这个时候只需要用单个exon的注释就可以, 从里面选一个;
                    // 如果也没有, 那么就是NA;
                    //单个exon的情况;
                    TempDist.clear();
                    for (const auto& every_gene_tx:SingleExonAnno) {
                        if (every_gene_tx.second[0] <= BE[1] && BE[0] <= every_gene_tx.second[1]) {
                            size_t pos = every_gene_tx.first.find('|');
                            if (pos != std::string::npos) {
                                // 根据 "|" 将字符串分割为两部分
                                first_part = every_gene_tx.first.substr(0, pos);
                            }
                            TempGene.push_back(first_part);
                            int a = abs(BE[1] - every_gene_tx.second[1]);
                            int b = abs(every_gene_tx.second[0] - BE[0]);
                            TempDist.push_back((a>b)?a:b);                            
                        }
                    }
                    if (TempDist.size() > 0) {
                        auto min_it = std::min_element(TempDist.begin(), TempDist.end());
                        first_part = TempGene[std::distance(TempDist.begin(), min_it)];
                    } else {
                        first_part = "NA";
                    }                    
                }

                FinalAnnotations.transcript2gene[itsname] = first_part;
                {
                    std::unique_lock<std::mutex> lock(updatedGtfMutex);
                    Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "transcript" << '\t' << BE[0] << '\t' << BE[1] << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\";" << '\n';
                    for (int i = 0; i < itssj.size()+1; i++){
                        if (i == 0) {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << BE[0] << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" <<  '\n';
                        } else if (i == itssj.size()) {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << BE[1] << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                        } else {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                        }
                    }
                }
                //输出追踪文件;
                if (Trace.is_open()){
                    for (const auto& EachRead:HighClusters[nameya]){
                        std::unique_lock<std::mutex> lock(traceMutex);
                        Trace << EachRead << '\t' << "novel_isoform" << '\t' << itsname << '\t' << first_part << '\t' << groupreadfiles[EachRead] << '\n'; 
                    }
                }
            }
        }
    }
    return FinalAnnotations;
}

OutputInformation Write_Detection_Transcript2gtf_MultiFiles(std::ofstream& Updated_Files, std::ofstream& Trace,
                                                DetectionResults& noderesults, DistanceInform& Disinform, 
                                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupannotations, 
                                                std::map<std::string, std::vector<std::array<int,2>>>& Annoexon, 
                                                std::unordered_map<std::size_t, std::array<int, 2>>& groupreadcoverage, 
                                                std::set<std::string>& groupallgene, 
                                                std::map<std::string, std::array<int,2>>& Annogenecovergae, 
                                                std::unordered_map<std::string, std::string>& GTF_Transcript_Strand, 
                                                std::unordered_map<size_t, std::string>& High_Strand, 
                                                std::string& chrname, std::string& group_size, 
                                                std::map<size_t, std::vector<std::string>>& HighClusters,
                                                int& FileNumber, std::unordered_map<std::string, std::string>& groupreadfiles,
                                                Solvent& filesolvent,
                                                std::unordered_map<std::string, std::array<int,2>>& SingleExonAnno,
                                                std::map<std::string, std::vector<std::string>>& gtf_gene2tx){
    OutputInformation FinalAnnotations;
    std::string itsname;
    std::vector<std::array<int,2>> itssj;
    std::vector<std::array<int,2>> itsexon;
    std::vector<std::string> TempGene;
    std::string TempGeneName;
    std::vector<int> TempDist;
    int novel_count = 0;
    std::size_t nameya;
    std::array<int,2> BE;
    std::string first_part;
    std::string second_part;

    for (int i = 0; i < FileNumber; i++) {
        FinalAnnotations.File_TranscriptNumber[std::to_string(i)] = {};
    }

    if (Updated_Files.is_open()) {
        // 第一次
        for (const auto aaa:noderesults.TrueNodeSet){
            if (aaa < Disinform.Index2Anno.size()){
                //有注释的部分;
                itsname = Disinform.Index2Anno[aaa];
                itssj = groupannotations[itsname];
                itsexon = Annoexon[itsname];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];

                for (int i = 0; i < FileNumber; i++) {
                    auto is1 = filesolvent.File_FSM[std::to_string(i)].find(itsname);
                    if (is1 != filesolvent.File_FSM[std::to_string(i)].end()) {
                        FinalAnnotations.File_TranscriptNumber[std::to_string(i)][itsname] = filesolvent.File_FSM[std::to_string(i)][itsname];
                    } else {
                        FinalAnnotations.File_TranscriptNumber[std::to_string(i)][itsname] = 0;
                    }
                }

                size_t pos = itsname.find('|');
                if (pos != std::string::npos) {
                    // 根据 "|" 将字符串分割为两部分
                    first_part = itsname.substr(0, pos);
                    second_part = itsname.substr(pos + 1);
                }
                {
                    std::unique_lock<std::mutex> lock(updatedGtfMutex);
                    Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t' << itsexon[0][0] << '\t' << itsexon[itsexon.size()-1][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\";" << '\n';
                    for (int i = 0; i < itsexon.size(); i++){
                        Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t' << itsexon[i][0] << '\t' << itsexon[i][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\"; exon_number \"" << i << "\";" << '\n';
                    }
                }
                FinalAnnotations.transcript2gene[itsname] = first_part;

            } else {
                // 这个是novel了;
                novel_count = novel_count + 1;
                //没有注释的部分;
                itsname = chrname + "-novel-" + group_size + "-" + std::to_string(aaa) + "-" + std::to_string(novel_count);
                Disinform.Index2novelname[aaa] = itsname;
                itssj = Disinform.Index2Unknown[aaa];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];
                nameya = Disinform.Index2hashname[aaa];
                BE = groupreadcoverage[nameya];
                
                for (int i = 0; i < FileNumber; i++) {
                    auto is1 = filesolvent.File_HighConClusters[std::to_string(i)].find(nameya);
                    if (is1 != filesolvent.File_HighConClusters[std::to_string(i)].end()) {
                        FinalAnnotations.File_TranscriptNumber[std::to_string(i)][itsname] = filesolvent.File_HighConClusters[std::to_string(i)][nameya].size();
                    } else {
                        FinalAnnotations.File_TranscriptNumber[std::to_string(i)][itsname] = 0;
                    }
                }
                
                TempGene.clear();
                // gene set里面已经有很多isoform; 当前isoform的sj就是itssj;
                if (groupallgene.size() != 0) {
                    TempGeneName = get_co_sj_gene(itssj, groupallgene, gtf_gene2tx, groupannotations);

                    if (TempGeneName.size() != 0) {
                        first_part = TempGeneName;
                    } else {
                        // 这个时候就是多个exon的注释和单个exon的注释, 从里面选一个;
                        for (const auto& every_gene:groupallgene) {
                            if (Annogenecovergae[every_gene][0] <= BE[1] && BE[0] <= Annogenecovergae[every_gene][1]) {
                                TempGene.push_back(every_gene);
                            }
                        }
                        //novel isoform没有候选的基因;
                        if (TempGene.size() == 0) {
                            first_part = "NA";
                            // novel isoform只有1个候选的基因;
                        } else if (TempGene.size() == 1) {
                            first_part = TempGene[0];
                            //novel isoform有2个以上候选的基因;
                        } else {
                            TempDist.clear();
                            for (const auto& each_gene:TempGene){
                                int a = abs(BE[1] - Annogenecovergae[each_gene][1]);
                                int b = abs(Annogenecovergae[each_gene][0] - BE[0]);
                                TempDist.push_back((a>b)?a:b);
                            }   
                            auto min_it = std::min_element(TempDist.begin(), TempDist.end());
                            first_part = TempGene[std::distance(TempDist.begin(), min_it)];                           
                        }
                    }

                } else {
 
                    TempDist.clear();
                    for (const auto& every_gene_tx:SingleExonAnno) {
                        if (every_gene_tx.second[0] <= BE[1] && BE[0] <= every_gene_tx.second[1]) {
                            size_t pos = every_gene_tx.first.find('|');
                            if (pos != std::string::npos) {
                                // 根据 "|" 将字符串分割为两部分
                                first_part = every_gene_tx.first.substr(0, pos);
                            }
                            TempGene.push_back(first_part);
                            int a = abs(BE[1] - every_gene_tx.second[1]);
                            int b = abs(every_gene_tx.second[0] - BE[0]);
                            TempDist.push_back((a>b)?a:b);                            
                        }
                    }
                    if (TempDist.size() > 0) {
                        auto min_it = std::min_element(TempDist.begin(), TempDist.end());
                        first_part = TempGene[std::distance(TempDist.begin(), min_it)];
                    } else {
                        first_part = "NA";
                    }                    
                }

                FinalAnnotations.transcript2gene[itsname] = first_part;
                // std::cout << first_part << std::endl;
                {
                    std::unique_lock<std::mutex> lock(updatedGtfMutex);
                    Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "transcript" << '\t' << BE[0] << '\t' << BE[1] << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\";" << '\n';
                    for (int i = 0; i < itssj.size()+1; i++){
                        if (i == 0) {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << BE[0] << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" <<  '\n';
                        } else if (i == itssj.size()) {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << BE[1] << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                        } else {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                        }
                    }
                }
                //输出追踪文件;
                if (Trace.is_open()){
                    for (const auto& EachRead:HighClusters[nameya]){
                        std::unique_lock<std::mutex> lock(traceMutex);
                        Trace << EachRead << '\t' << "novel_isoform" << '\t' << itsname << '\t' << first_part << '\t' << groupreadfiles[EachRead] << '\n'; 
                    }
                }
            }
        }
    }
    return FinalAnnotations;
}


struct FileFalseNode
{
    std::set<int> ThisFlaseNode;
    std::unordered_map<int, std::vector<std::string>> FasleNode2Count;
};
FileFalseNode get_File_False_node (std::set<int>& FalseNodeSet, 
                          std::unordered_map<int, std::size_t>& Node2hashname, 
                          std::map<std::string, std::map<size_t, std::vector<std::string>>>& File_HighClusters, 
                          int & FileNNumber) {
    //这个Barcode中所有错误的节点;
    FileFalseNode ThisFile_FalseNumber;
    std::map<size_t, std::vector<std::string>> this_file_HC = File_HighClusters[std::to_string(FileNNumber)];
    for (const auto& node:FalseNodeSet) {
        //每个错误节点不一定在这个File的HighCluster中;
        auto which = this_file_HC.find(Node2hashname[node]);
        if (which != this_file_HC.end()) {
            ThisFile_FalseNumber.ThisFlaseNode.insert(node);
            ThisFile_FalseNumber.FasleNode2Count[node] = this_file_HC[Node2hashname[node]];
        }
    }
    return ThisFile_FalseNumber;
}


int whether_isoform_part_is_not(std::vector<std::array<int,2>> AnnoSJ, std::vector<std::array<int,2>> ISMSJ){
    int ISM_Flag = 0;
    auto it1 = std::find(AnnoSJ.begin(), AnnoSJ.end(), ISMSJ[0]);
    if (it1 != AnnoSJ.end()){
        //计算出it1的索引;
        int it1_index = it1 - AnnoSJ.begin();
        // std::cout << it1_index << std::endl;
        auto it2 = std::find(AnnoSJ.begin(), AnnoSJ.end(), ISMSJ[ISMSJ.size()-1]);
        if (it2 != AnnoSJ.end()){
            //计算出it2的索引;
            int it2_index = it2 - AnnoSJ.begin();
            // std::cout << it2_index << std::endl;
            //如果两个索引之间相差的大小等于长度减1;
            if (it2_index - it1_index == ISMSJ.size() - 1){
                ISM_Flag = 1;
            }
        }
    } 
    return ISM_Flag;
}

std::string concatenateSet(const std::set<std::string>& stringSet) {
    if (stringSet.empty()) {
        return ""; // 如果集合为空，返回空字符串
    }
    std::ostringstream oss;
    std::string oss_string;
    for (const auto& str : stringSet) {
        oss << str << ',';
    }
    oss_string = oss.str();
    oss_string.pop_back();
    return oss_string;
}

struct IndicateFire
{
    Eigen::MatrixXd Indicate_Matrix;
    Eigen::VectorXd Cluster_Number;
    std::vector<std::string> Order_Transcript_Name_Vector;
};


IndicateFire Quantification_initialization(std::map<std::size_t, std::vector<std::string>>& groupISM, 
                                           OutputInformation& FinallyAnnotations, 
                                           std::set<int>& truenodeset,
                                           std::set<int>& falsenodeset, 
                                           DistanceInform& Disinform, 
                                           std::string& group_size, 
                                           std::map<size_t, std::vector<std::string>>& HighClusters, 
                                           std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs, 
                                           std::unordered_map<std::string, std::string>& groupreadfiles,
                                           std::ofstream& Trace,
                                           std::string& chrname) {
    // std::mutex& TraceMutex;
    IndicateFire OutputResults;
    //生成的变量, 识别的转录本按照顺序;
    std::vector<std::string> Order_Transcript_Name;
    //生成的变量, 最终的ISM指示矩阵;
    Eigen::MatrixXd known_ISM_matrix;
    std::vector<std::array<int,2>> AreadSJs;
    
    //建立索引和判断条件;
    int RowIndex = -1;
    int ColIndex = -1;
    int flag = 0;


    for (const auto& eachTransctipt:FinallyAnnotations.Transcript_Annotations){
        Order_Transcript_Name.push_back(eachTransctipt.first);
    }
    Eigen::RowVectorXd EachClusterRowVector(Order_Transcript_Name.size());
    //存储一个列向量, 每个类的个数, 方便后续做乘法;
    Eigen::VectorXd known_ISM_cluster_numbers;

    //ISM的追踪文件;
    std::vector<std::string> ism_gtf_name;
    std::set<std::string> first_part_set;
    std::set<std::string> second_part_set;
    std::string first_part;
    std::string second_part;
    std::vector<int> ColIndexVec;
    double max_ratio = 0;
    double this_ratio = 0;

    //这一段得到ISM的指示矩阵;
    if (groupISM.size() != 0){
        RowIndex = -1;
        for (const auto& eachISM:groupISM){
            //对生成的每个行向量清空; 每个cluster的行向量;
            AreadSJs = groupreadsjs[eachISM.second[0]];
            EachClusterRowVector.setZero();
            ColIndex = -1;
            ism_gtf_name.clear();
            first_part_set.clear();
            second_part_set.clear();
            ColIndexVec.clear();   
            double max_ratio = 0;
            double this_ratio = 0;
            //对每个注释遍历;
            for (const auto& eachTransctipt:Order_Transcript_Name){
                ColIndex = ColIndex + 1;
                //如果ISM的个数大于注释, 则不可能是1;
                if (AreadSJs.size() < FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size()){
                    if ((AreadSJs[0][0] >= FinallyAnnotations.Transcript_Annotations[eachTransctipt].first[0][0]) && (AreadSJs[AreadSJs.size()-1][1] <= FinallyAnnotations.Transcript_Annotations[eachTransctipt].first[(FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size()-1)][1])){
                        flag = whether_isoform_part_is_not(FinallyAnnotations.Transcript_Annotations[eachTransctipt].first, AreadSJs);
                        if (flag == 1){
                            this_ratio = static_cast<double>(AreadSJs.size())/FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size();
                            if (max_ratio < this_ratio) {
                                max_ratio = this_ratio;
                            }
                            EachClusterRowVector(ColIndex) = 1;
                            ism_gtf_name.push_back(eachTransctipt);
                            ColIndexVec.push_back(ColIndex);
                        } else {
                            EachClusterRowVector(ColIndex) = 0;
                        }
                    } else{
                        EachClusterRowVector(ColIndex) = 0;
                    }
                } else{
                    EachClusterRowVector(ColIndex) = 0;
                }
            }

            // 在这里加一步过滤最好了, ISM如果有一个特别近的, 还有一个特别远的, 特别远的就去掉吧;
            if (ColIndexVec.size() > 1) {
                if (max_ratio >= 0.5) {
                    for (const auto& ColIndex:ColIndexVec) {
                        std::string TransName = Order_Transcript_Name[ColIndex];
                        this_ratio = static_cast<double>(AreadSJs.size())/FinallyAnnotations.Transcript_Annotations[TransName].first.size();
                        if (this_ratio < max_ratio) {
                            EachClusterRowVector[ColIndex] = 0;
                            ism_gtf_name.erase(std::remove(ism_gtf_name.begin(), ism_gtf_name.end(), TransName), ism_gtf_name.end());
                        }
                    }
                }
            } 

            if (ism_gtf_name.size() != 0) {
                //只有这一行不全为0, 才能放入矩阵;
                RowIndex = RowIndex + 1;
                //调整矩阵大小以容纳第一个行向量;
                known_ISM_matrix.conservativeResize(RowIndex+1, Order_Transcript_Name.size());
                known_ISM_matrix.row(RowIndex) = EachClusterRowVector;
                //存储一个列向量, 每个类的个数, 方便后续做乘法;
                known_ISM_cluster_numbers.conservativeResize(RowIndex+1, 1);
                known_ISM_cluster_numbers(RowIndex, 0) = eachISM.second.size();

                //开始输出ISM的追踪文件;
                for (const auto& EachTranscript:ism_gtf_name){
                    size_t pos = EachTranscript.find('|');
                    if (pos != std::string::npos) {
                        first_part_set.insert(EachTranscript.substr(0, pos));
                        second_part_set.insert(EachTranscript.substr(pos + 1));
                    }
                }
                first_part = concatenateSet(first_part_set);
                second_part = concatenateSet(second_part_set);
                if (Trace.is_open()) {
                    for (const auto& EachRead:eachISM.second){
                        std::unique_lock<std::mutex> lock(traceMutex);
                        Trace << EachRead << '\t' << "ISM" << '\t' << second_part << '\t' << first_part << '\t' << groupreadfiles[EachRead] << '\n'; 
                    }
                }
            }
        }
    }

    std::set<int> MergeFalseSet(falsenodeset.begin(), falsenodeset.end());

    //初始化生成的指示矩阵;
    Eigen::MatrixXd False_novel_candidate_matrix(MergeFalseSet.size(), Order_Transcript_Name.size());
    //每个错误节点的所有邻居集合;
    std::set<int> FalseNode_NeighborSet;
    std::set<int> FalseNode_TrueSet;
    //对应正确节点的名称;
    std::string itsname;
    //要查找的子字符串;
    std::string subStr;
    //构造一个行向量;
    Eigen::RowVectorXd FalseNode_RowVector(Order_Transcript_Name.size());
    RowIndex = -1;
    int countC = 0;
    
    if (MergeFalseSet.size() != 0){
        //构造识别出的所有正确节点;
        std::set<int> MergedTrueSet(truenodeset.begin(), truenodeset.end());
        
        //遍历所有识别为错误的节点;
        //找到每个错误节点的邻居节点 --- FalseNode_NeighborSet;
        for (const auto& FalseNode:MergeFalseSet){    
            RowIndex = RowIndex + 1;        
            //找到与他相邻的节点;
            FalseNode_NeighborSet.clear();
            FalseNode_TrueSet.clear();            
            FalseNode_NeighborSet = Disinform.NodeDistanceSet[FalseNode];

            std::set_intersection(MergedTrueSet.begin(), MergedTrueSet.end(), FalseNode_NeighborSet.begin(), FalseNode_NeighborSet.end(),
            std::inserter(FalseNode_TrueSet, FalseNode_TrueSet.begin()));

            FalseNode_RowVector.setZero();
            ism_gtf_name.clear();

            if (FalseNode_TrueSet.size() != 0) {
                //对所有的正确节点遍历;
                for(const auto& EachTNode:FalseNode_TrueSet){
                    if (EachTNode < Disinform.Index2Anno.size()){
                        //说明这个邻居是有注释的转录本;
                        itsname = Disinform.Index2Anno[EachTNode];
                        //接下来找itsname在Order_Transcript_Name中的索引;
                        auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                        FalseNode_RowVector[itwhere] = 1;
                        ism_gtf_name.push_back(itsname);

                    } else{
                        //这个邻居是novel的转录本;
                        itsname = Disinform.Index2novelname[EachTNode];
                        auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                        FalseNode_RowVector[itwhere] = 1;
                        ism_gtf_name.push_back(itsname);
                    }
                }
                //将这一行加入到矩阵中;
                False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
            } 
            else
            {   
                countC = 0;
                //需要找邻居的邻居;
                for(const auto& neinode:FalseNode_NeighborSet){
                    if (Disinform.NodeDistanceSet[neinode].size() != 0){
                        for (const auto& node:Disinform.NodeDistanceSet[neinode]){
                            if (node < Disinform.Index2Anno.size()){
                                //说明这个邻居是有注释的转录本;
                                itsname = Disinform.Index2Anno[node];
                                //接下来找itsname在Order_Transcript_Name中的索引;
                                auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                                FalseNode_RowVector[itwhere] = 1;
                                ism_gtf_name.push_back(itsname);
                                countC++;
                            } else {
                                if (truenodeset.find(node) != truenodeset.end()) {
                                    // 这个邻居是novel的转录本;
                                    itsname = Disinform.Index2novelname[node];
                                    auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                                    FalseNode_RowVector[itwhere] = 1;
                                    ism_gtf_name.push_back(itsname);
                                    countC++;
                                }
                            }
                        }
                    }
                }
                if (countC != 0){
                    False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
                } else {
                    //这种情况下就是邻居的邻居也没有一个好的, 那就找距离最小值了;
                    std::vector<int> houX = {};
                    int distance = 0;
                    for (const auto& isofrorm:Order_Transcript_Name){
                        distance = IntervalMerge(FinallyAnnotations.Transcript_Annotations[isofrorm].first, Disinform.Index2Unknown[FalseNode]) - IntervalIntersection(FinallyAnnotations.Transcript_Annotations[isofrorm].first, Disinform.Index2Unknown[FalseNode]);
                        houX.push_back(distance);
                    }
                    // 使用 std::min_element 找到最小元素的迭代器
                    auto min_it = std::min_element(houX.begin(), houX.end());

                    // 计算最小元素的索引
                    int min_index = std::distance(houX.begin(), min_it);
                    FalseNode_RowVector[min_index] = 1;
                    False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
                    ism_gtf_name.push_back(Order_Transcript_Name[min_index]);
                    // std::cout << "*** Neighbors neighbors are wrong. That doesn't happen!!! ***" << std::endl;
                }
            }
            first_part_set.clear();
            second_part_set.clear();
            //开始输出HighCon的追踪文件;
            for (const auto& EachTranscript:ism_gtf_name){
                size_t pos = EachTranscript.find('|');
                if (pos != std::string::npos) {
                    first_part_set.insert(EachTranscript.substr(0, pos));
                    second_part_set.insert(EachTranscript.substr(pos + 1));
                } else {
                    first_part_set.insert(FinallyAnnotations.transcript2gene[EachTranscript]);
                    second_part_set.insert(EachTranscript);
                }
            }
            first_part = concatenateSet(first_part_set);
            second_part = concatenateSet(second_part_set);
            if (Trace.is_open()){
                for(const auto& EachRead:HighClusters[Disinform.Index2hashname[FalseNode]]){
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << EachRead << '\t' << "approximate" << '\t' << second_part << '\t' << first_part << '\t' << groupreadfiles[EachRead] << '\n';
                }
            }            
        }
    }

    Eigen::VectorXd False_novel_candidate_cluster_numbers(MergeFalseSet.size());
    RowIndex = -1;
    for (const auto& FalseNode:MergeFalseSet){
        RowIndex = RowIndex + 1;
        False_novel_candidate_cluster_numbers[RowIndex] = Disinform.Index2Count[FalseNode];
    }


    //最后输出结果矩阵和结果量;
    if (known_ISM_matrix.rows() != 0){
        if (False_novel_candidate_matrix.rows() != 0){
            OutputResults.Indicate_Matrix.resize(known_ISM_matrix.rows()+False_novel_candidate_matrix.rows(), known_ISM_matrix.cols());
            OutputResults.Cluster_Number.resize(known_ISM_matrix.rows()+False_novel_candidate_matrix.rows(), known_ISM_cluster_numbers.cols());
            OutputResults.Indicate_Matrix << known_ISM_matrix,
                                             False_novel_candidate_matrix;
            OutputResults.Cluster_Number << known_ISM_cluster_numbers,
                                            False_novel_candidate_cluster_numbers;
        }
        else
        {
            OutputResults.Indicate_Matrix = known_ISM_matrix;
            OutputResults.Cluster_Number = known_ISM_cluster_numbers;
        }
    } else{
        if (False_novel_candidate_matrix.rows() != 0){
            OutputResults.Indicate_Matrix = False_novel_candidate_matrix;
            OutputResults.Cluster_Number = False_novel_candidate_cluster_numbers;
        }
    }
    OutputResults.Order_Transcript_Name_Vector = Order_Transcript_Name;
    return OutputResults;
}



IndicateFire Quantification_initialization_MultiFiles (std::map<std::size_t, std::vector<std::string>>& groupISM, 
                                                OutputInformation& FinallyAnnotations, 
                                                std::set<int>& truenodeset,
                                                std::set<int>& falsenodeset,
                                                std::unordered_map<int, std::vector<std::string>>& faslenode2count, 
                                                DistanceInform& Disinform, 
                                                std::string& group_size, 
                                                std::map<size_t, std::vector<std::string>>& HighClusters, 
                                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs,
                                                std::unordered_map<std::string, std::string>& groupreadfiles,
                                                std::ofstream& Trace) {
    // std::mutex& TraceMutex;
    IndicateFire OutputResults;
    //生成的变量, 识别的转录本按照顺序;
    std::vector<std::string> Order_Transcript_Name;
    //生成的变量, 最终的ISM指示矩阵;
    Eigen::MatrixXd known_ISM_matrix;
    std::vector<std::array<int,2>> AreadSJs;
    
    //建立索引和判断条件;
    int RowIndex = -1;
    int ColIndex = -1;
    int flag = 0;

    for (const auto& eachTransctipt:FinallyAnnotations.Transcript_Annotations) {
        Order_Transcript_Name.push_back(eachTransctipt.first);
    }
    Eigen::RowVectorXd EachClusterRowVector(Order_Transcript_Name.size());
    //存储一个列向量, 每个类的个数, 方便后续做乘法;
    Eigen::VectorXd known_ISM_cluster_numbers;

    //ISM的追踪文件;
    std::vector<std::string> ism_gtf_name;
    std::set<std::string> first_part_set;
    std::set<std::string> second_part_set;
    std::string first_part;
    std::string second_part;
    std::vector<int> ColIndexVec;
    double max_ratio = 0;
    double this_ratio = 0;

    if (groupISM.size() != 0){
        RowIndex = -1;
        for (const auto& eachISM:groupISM){
            AreadSJs = groupreadsjs[eachISM.second[0]];
            EachClusterRowVector.setZero();
            ColIndex = -1;
            ism_gtf_name.clear();
            first_part_set.clear();
            second_part_set.clear();            
            ColIndexVec.clear();   
            double max_ratio = 0;
            double this_ratio = 0;
            //对每个注释遍历;
            for (const auto& eachTransctipt:Order_Transcript_Name){
                ColIndex = ColIndex + 1;
                //如果ISM的个数大于注释, 则不可能是1;
                if (AreadSJs.size() < FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size()){
                    if ((AreadSJs[0][0] >= FinallyAnnotations.Transcript_Annotations[eachTransctipt].first[0][0]) && (AreadSJs[AreadSJs.size()-1][1] <= FinallyAnnotations.Transcript_Annotations[eachTransctipt].first[(FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size()-1)][1])){
                        flag = whether_isoform_part_is_not(FinallyAnnotations.Transcript_Annotations[eachTransctipt].first, AreadSJs);
                        if (flag == 1){
                            this_ratio = static_cast<double>(AreadSJs.size())/FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size();
                            if (max_ratio < this_ratio) {
                                max_ratio = this_ratio;
                            }
                            EachClusterRowVector(ColIndex) = 1;
                            ism_gtf_name.push_back(eachTransctipt);
                            ColIndexVec.push_back(ColIndex);
                        } else {
                            EachClusterRowVector(ColIndex) = 0;
                        }
                    } else{
                        EachClusterRowVector(ColIndex) = 0;
                    }
                } else {
                    EachClusterRowVector(ColIndex) = 0;
                }
            }
            // 在这里加一步过滤最好了, ISM如果有一个特别近的, 还有一个特别远的, 特别远的就去掉吧;
            if (ColIndexVec.size() > 1) {
                if (max_ratio >= 0.5) {
                    for (const auto& ColIndex:ColIndexVec) {
                        std::string TransName = Order_Transcript_Name[ColIndex];
                        this_ratio = static_cast<double>(AreadSJs.size())/FinallyAnnotations.Transcript_Annotations[TransName].first.size();
                        if (this_ratio < max_ratio) {
                            EachClusterRowVector[ColIndex] = 0;
                            ism_gtf_name.erase(std::remove(ism_gtf_name.begin(), ism_gtf_name.end(), TransName), ism_gtf_name.end());
                            // std::cout << " 删除 " << TransName << std::endl;
                        }
                    }
                }
            } 
            if (ism_gtf_name.size() != 0){
                //只有这一行不全为0, 才能放入矩阵;
                RowIndex = RowIndex + 1;
                //调整矩阵大小以容纳第一个行向量;
                known_ISM_matrix.conservativeResize(RowIndex+1, Order_Transcript_Name.size());
                known_ISM_matrix.row(RowIndex) = EachClusterRowVector;
                //存储一个列向量, 每个类的个数, 方便后续做乘法;
                known_ISM_cluster_numbers.conservativeResize(RowIndex+1, 1);
                known_ISM_cluster_numbers(RowIndex, 0) = eachISM.second.size();

                //开始输出ISM的追踪文件;
                for (const auto& EachTranscript:ism_gtf_name){
                    size_t pos = EachTranscript.find('|');
                    if (pos != std::string::npos) {
                        first_part_set.insert(EachTranscript.substr(0, pos));
                        second_part_set.insert(EachTranscript.substr(pos + 1));
                    }
                }
                first_part = concatenateSet(first_part_set);
                second_part = concatenateSet(second_part_set);
                if (Trace.is_open()){
                    for (const auto& EachRead:eachISM.second){
                        std::unique_lock<std::mutex> lock(traceMutex);
                        Trace << EachRead << '\t' << "ISM" << '\t' << second_part << '\t' << first_part << '\t' << groupreadfiles[EachRead] << '\n'; 
                    }
                }
            } 
        }
    }
    std::set<int> MergeFalseSet(falsenodeset.begin(), falsenodeset.end());

    //初始化生成的指示矩阵;
    Eigen::MatrixXd False_novel_candidate_matrix(MergeFalseSet.size(), Order_Transcript_Name.size());
    //每个错误节点的所有邻居集合;
    std::set<int> FalseNode_NeighborSet;
    std::set<int> FalseNode_TrueSet;
    //对应正确节点的名称;
    std::string itsname;
    //要查找的子字符串;
    std::string subStr;
    //构造一个行向量;
    Eigen::RowVectorXd FalseNode_RowVector(Order_Transcript_Name.size());
    RowIndex = -1;
    int countC = 0;
    
    if (MergeFalseSet.size() != 0) {
        //构造识别出的所有正确节点;
        std::set<int> MergedTrueSet(truenodeset.begin(), truenodeset.end());
        
        //遍历所有识别为错误的节点;
        //找到每个错误节点的邻居节点 --- FalseNode_NeighborSet;
        for (const auto& FalseNode:MergeFalseSet){    
            RowIndex = RowIndex + 1;        
            //找到与他相邻的节点;
            FalseNode_NeighborSet.clear();
            FalseNode_TrueSet.clear();            
            FalseNode_NeighborSet = Disinform.NodeDistanceSet[FalseNode];

            std::set_intersection(MergedTrueSet.begin(), MergedTrueSet.end(), FalseNode_NeighborSet.begin(), FalseNode_NeighborSet.end(),
            std::inserter(FalseNode_TrueSet, FalseNode_TrueSet.begin()));

            FalseNode_RowVector.setZero();
            ism_gtf_name.clear();

            if (FalseNode_TrueSet.size() != 0) {
                //对所有的正确节点遍历;
                for(const auto& EachTNode:FalseNode_TrueSet){
                    if (EachTNode < Disinform.Index2Anno.size()){
                        //说明这个邻居是有注释的转录本;
                        itsname = Disinform.Index2Anno[EachTNode];
                        //接下来找itsname在Order_Transcript_Name中的索引;
                        auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                        FalseNode_RowVector[itwhere] = 1;
                        ism_gtf_name.push_back(itsname);

                    } else{
                        //这个邻居是novel的转录本;
                        itsname = Disinform.Index2novelname[EachTNode];
                        auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                        FalseNode_RowVector[itwhere] = 1;
                        ism_gtf_name.push_back(itsname);
                    }
                }
                //将这一行加入到矩阵中;
                False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
            } 
            else
            {   
                countC = 0;
                //需要找邻居的邻居;
                for(const auto& neinode:FalseNode_NeighborSet){
                    if (Disinform.NodeDistanceSet[neinode].size() != 0){
                        for (const auto& node:Disinform.NodeDistanceSet[neinode]){
                            if (node < Disinform.Index2Anno.size()){
                                //说明这个邻居是有注释的转录本;
                                itsname = Disinform.Index2Anno[node];
                                //接下来找itsname在Order_Transcript_Name中的索引;
                                auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                                FalseNode_RowVector[itwhere] = 1;
                                ism_gtf_name.push_back(itsname);
                                countC++;
                            } else {
                                if (truenodeset.find(node) != truenodeset.end()) {
                                    // 这个邻居是novel的转录本;
                                    itsname = Disinform.Index2novelname[node];
                                    auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                                    FalseNode_RowVector[itwhere] = 1;
                                    ism_gtf_name.push_back(itsname);
                                    countC++;
                                }
                            }
                        }
                    }
                }
                if (countC != 0){
                    False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
                } else {
                    //这种情况下就是邻居的邻居也没有一个好的, 那就找距离最小值了;
                    std::vector<int> houX = {};
                    int distance = 0;
                    for (const auto& isofrorm:Order_Transcript_Name){
                        distance = IntervalMerge(FinallyAnnotations.Transcript_Annotations[isofrorm].first, Disinform.Index2Unknown[FalseNode]) - IntervalIntersection(FinallyAnnotations.Transcript_Annotations[isofrorm].first, Disinform.Index2Unknown[FalseNode]);
                        houX.push_back(distance);
                    }
                    // 使用 std::min_element 找到最小元素的迭代器
                    auto min_it = std::min_element(houX.begin(), houX.end());

                    // 计算最小元素的索引
                    int min_index = std::distance(houX.begin(), min_it);
                    FalseNode_RowVector[min_index] = 1;
                    False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
                    ism_gtf_name.push_back(Order_Transcript_Name[min_index]);
                    // std::cout << "*** Neighbors neighbors are wrong. That doesn't happen!!! ***" << std::endl;
                }
            }
            first_part_set.clear();
            second_part_set.clear();
            //开始输出HighCon的追踪文件;
            for (const auto& EachTranscript:ism_gtf_name){
                size_t pos = EachTranscript.find('|');
                if (pos != std::string::npos) {
                    first_part_set.insert(EachTranscript.substr(0, pos));
                    second_part_set.insert(EachTranscript.substr(pos + 1));
                } else {
                    first_part_set.insert(FinallyAnnotations.transcript2gene[EachTranscript]);
                    second_part_set.insert(EachTranscript);
                }
            }
            first_part = concatenateSet(first_part_set);
            second_part = concatenateSet(second_part_set);
            if (Trace.is_open()){
                for(const auto& EachRead:HighClusters[Disinform.Index2hashname[FalseNode]]){
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << EachRead << '\t' << "approximate" << '\t' << second_part << '\t' << first_part << '\t' << groupreadfiles[EachRead] << '\n';
                }
            }            
        }
    }
    //存储一个列向量, 每个类的个数, 方便后续做乘法;
    Eigen::VectorXd False_novel_candidate_cluster_numbers(MergeFalseSet.size());
    RowIndex = -1;
    for (const auto& FalseNode:MergeFalseSet){
        RowIndex = RowIndex + 1;
        False_novel_candidate_cluster_numbers[RowIndex] = faslenode2count[FalseNode].size();
    }

    //最后输出结果矩阵和结果量;
    if (known_ISM_matrix.rows() != 0){
        if (False_novel_candidate_matrix.rows() != 0){
            OutputResults.Indicate_Matrix.resize(known_ISM_matrix.rows()+False_novel_candidate_matrix.rows(), known_ISM_matrix.cols());
            OutputResults.Cluster_Number.resize(known_ISM_matrix.rows()+False_novel_candidate_matrix.rows(), known_ISM_cluster_numbers.cols());
            OutputResults.Indicate_Matrix << known_ISM_matrix,
                                             False_novel_candidate_matrix;
            OutputResults.Cluster_Number << known_ISM_cluster_numbers,
                                            False_novel_candidate_cluster_numbers;
        }
        else
        {
            OutputResults.Indicate_Matrix = known_ISM_matrix;
            OutputResults.Cluster_Number = known_ISM_cluster_numbers;
        }
    } else{
        if (False_novel_candidate_matrix.rows() != 0){
            OutputResults.Indicate_Matrix = False_novel_candidate_matrix;
            OutputResults.Cluster_Number = False_novel_candidate_cluster_numbers;
        }
    }
    OutputResults.Order_Transcript_Name_Vector = Order_Transcript_Name;
    return OutputResults;
}


std::map<std::string, double> EM_Alg (OutputInformation& FinallyAnnotations, IndicateFire& Indicate_Number, 
             std::ofstream& isoformCountPath) {
    // 输出;
    std::map<std::string, double> geneCounts;
    //直接输出量化结果;
    if (Indicate_Number.Indicate_Matrix.rows() != 0 && Indicate_Number.Indicate_Matrix.cols() != 0){
        std::string its_name;
        // 初始化概率向量
        Eigen::VectorXd P_Col_init0 = Eigen::VectorXd::Constant(
            Indicate_Number.Order_Transcript_Name_Vector.size(),
            1.0 / Indicate_Number.Order_Transcript_Name_Vector.size() 
        );
        // 初始化;
        Eigen::MatrixXd Z(Indicate_Number.Indicate_Matrix.rows(), Indicate_Number.Indicate_Matrix.cols());
        Eigen::VectorXd P1(P_Col_init0.size());
        Eigen::MatrixXd AnnoN;
        int CountCyc = 1; // 初始化为 1，与原代码对齐
        double sum_abs_diff;
        do {

            Z = Indicate_Number.Indicate_Matrix.array().rowwise() * P_Col_init0.transpose().array();
            Z.array().colwise() /= (Indicate_Number.Indicate_Matrix * P_Col_init0).array();

            AnnoN = ((Indicate_Number.Cluster_Number.transpose()) * Z).transpose().reshaped(Indicate_Number.Indicate_Matrix.cols(),1);

            P1 = AnnoN / AnnoN.sum();

            sum_abs_diff = (P1 - P_Col_init0).cwiseAbs().sum();

            P_Col_init0 = P1;
            CountCyc = CountCyc + 1;
        } while ((sum_abs_diff > 5e-2) || (CountCyc <= 10));
    
        for (int i = 0; i < AnnoN.rows(); i++){
            its_name = Indicate_Number.Order_Transcript_Name_Vector[i];
            FinallyAnnotations.Transcript_Annotations[its_name].second = FinallyAnnotations.Transcript_Annotations[its_name].second + AnnoN(i,0);
        }
    }

    std::string geneName;
    std::string second_part;
    std::string first_part;

    if (isoformCountPath.is_open()){
        // 这里是多个exon的;
        for (const auto& eachAnno:FinallyAnnotations.Transcript_Annotations){
            size_t pos = eachAnno.first.find('|');
            if (pos != std::string::npos) {
                // 根据 "|" 将字符串分割为两部分
                second_part = eachAnno.first.substr(pos + 1);
            } else {
                second_part = eachAnno.first;
            }
            geneName = FinallyAnnotations.transcript2gene[eachAnno.first];
            if (eachAnno.second.second > 0){
                std::unique_lock<std::mutex> lock(isoformCountMutex);
                isoformCountPath << second_part << '\t' << geneName << '\t' << int(eachAnno.second.second) << '\n';
            }
            if (geneName != "NA"){
                auto it = geneCounts.find(FinallyAnnotations.transcript2gene[eachAnno.first]);
                if (it != geneCounts.end()){
                    geneCounts[geneName] = geneCounts[geneName] + eachAnno.second.second;
                } else {
                    geneCounts[geneName] = eachAnno.second.second;
                }
            }
        }
    }
    return geneCounts;
}   

// 最终EM算法的步骤;
void EM_Alg_MultiFiles (OutputInformation& FinallyAnnotations, 
                        IndicateFire& Indicate_Number, 
                        int& FileNumber) {
    //直接输出量化结果;
    if (Indicate_Number.Indicate_Matrix.rows() != 0 && Indicate_Number.Indicate_Matrix.cols() != 0) {
        // 初始化概率向量
        Eigen::VectorXd P_Col_init0 = Eigen::VectorXd::Constant(
            Indicate_Number.Order_Transcript_Name_Vector.size(),
            1.0 / Indicate_Number.Order_Transcript_Name_Vector.size() 
        );
        // 初始化;
        Eigen::MatrixXd Z(Indicate_Number.Indicate_Matrix.rows(), Indicate_Number.Indicate_Matrix.cols());
        Eigen::VectorXd P1(P_Col_init0.size());
        Eigen::MatrixXd AnnoN;
        int CountCyc = 1; 
        double sum_abs_diff;

        do {

            Z = Indicate_Number.Indicate_Matrix.array().rowwise() * P_Col_init0.transpose().array();
            Z.array().colwise() /= (Indicate_Number.Indicate_Matrix * P_Col_init0).array();

            AnnoN = ((Indicate_Number.Cluster_Number.transpose()) * Z).transpose().reshaped(Indicate_Number.Indicate_Matrix.cols(),1);

            P1 = AnnoN / AnnoN.sum();

            sum_abs_diff = (P1 - P_Col_init0).cwiseAbs().sum();

            P_Col_init0 = P1;
            CountCyc = CountCyc + 1;
        } while ((sum_abs_diff > 5e-2) || (CountCyc <= 10));
    

        std::string its_name;
        for (int i = 0; i < AnnoN.rows(); i++){
            its_name = Indicate_Number.Order_Transcript_Name_Vector[i];
            FinallyAnnotations.File_TranscriptNumber[std::to_string(FileNumber)][its_name] = FinallyAnnotations.File_TranscriptNumber[std::to_string(FileNumber)][its_name] + AnnoN(i,0);
        }
    }
}   

std::map<std::string, std::map<std::string, double>> get_allfile_genecounts(
                                            OutputInformation& FinallyAnnotations, 
                                            int& AllFileNumber, 
                                            std::ofstream& isoformCountPath){
    std::map<std::string, std::map<std::string, double>> geneCounts;
    std::string geneName;
    std::string second_part;
    //将结果输出到丰度文件; 这里是多个exon的;
    if (isoformCountPath.is_open()) {
        for (const auto& eachAnno:FinallyAnnotations.File_TranscriptNumber["0"]) {
            // 一个一个注释和novel;
            size_t pos = eachAnno.first.find('|');
            if (pos != std::string::npos) {
                // 根据 "|" 将字符串分割为两部分
                second_part = eachAnno.first.substr(pos + 1);
            } else {
                second_part = eachAnno.first;
            }
            geneName = FinallyAnnotations.transcript2gene[eachAnno.first];
            {   
                std::unique_lock<std::mutex> lock(isoformCountMutex);
                isoformCountPath << second_part << '\t' << geneName;
                for (int i = 0; i < AllFileNumber; i++) {
                    std::string filenumber = std::to_string(i);
                    isoformCountPath << '\t' << int(FinallyAnnotations.File_TranscriptNumber[filenumber][eachAnno.first]);
                }            
                isoformCountPath << '\n';
            } 
            if (geneName != "NA") {
                auto it = geneCounts.find(geneName);
                if (it != geneCounts.end()) {
                    for (int i = 0; i < AllFileNumber; i++) {
                        std::string filenumber = std::to_string(i);
                        geneCounts[geneName][filenumber] = geneCounts[geneName][filenumber] + FinallyAnnotations.File_TranscriptNumber[filenumber][eachAnno.first];
                    }
                } else {
                    geneCounts[geneName] = {};
                    for (int i = 0; i < AllFileNumber; i++) {
                        std::string filenumber = std::to_string(i);
                        geneCounts[geneName][filenumber] = FinallyAnnotations.File_TranscriptNumber[filenumber][eachAnno.first];
                    }                    
                }
            }           

        }
    }
    return geneCounts;
}

std::map<std::string, std::map<std::string, double>> DetectQuant(GroupAnnotation& groupanno, 
                 SpliceChainClass& splicechainclass, 
                 GroupInformation& groupinform, 
                 const int& groupdistance,
                 GTF& gtfexon, GTFsj& gtfsjs,
                 std::ofstream& gtfFilePath, std::ofstream& isoformFilePath,
                 std::ofstream& traceFilePath, int& FileNo,
                 std::unordered_map<std::string, std::array<int,2>>& singleexonanno) {
    // std::mutex& genemutex, std::mutex& isoformmutex, std::mutex& gtfmutex, std::mutex& tracemutex;
    // 识别和定量;
    std::string chrchr = groupinform.chrName;
    std::string groupnumber = groupinform.GroupIndex;
    DistanceInform DMatrix_GraphNode = get_distance_matrix(groupanno.Group_Annotations, 
                                                           splicechainclass.HighConClusters, 
                                                           groupinform.GroupReadSjs, 
                                                           groupdistance, splicechainclass.FSM);
    // 最终输出;
    std::map<std::string, std::map<std::string, double>> AllFile_GeneCounts;
    if (DMatrix_GraphNode.NodeDistanceSet.size() != 0) {
        // 寻找最大团;
        std::vector<std::vector<int>> CliquesVector = Find_Maximal_Cliques(DMatrix_GraphNode.NodeDistanceSet);

        // 识别;
        DetectionResults NodeResults = Transcript_Detection(CliquesVector, 
                                                            DMatrix_GraphNode.Index2Anno, 
                                                            DMatrix_GraphNode.Index2Unknown, 
                                                            DMatrix_GraphNode.Index2Count);
        // 量化前拆分;
        if (FileNo > 1) {
            // 多于一个sam文件, 就需要拆分;
            Solvent SpliceChainSolvent = get_Solvent_FsmIsmHigh(splicechainclass, FileNo, groupinform);
            // 写到更新的注释中;
            OutputInformation Finally_Annotations = Write_Detection_Transcript2gtf_MultiFiles(gtfFilePath, traceFilePath, 
                                                    NodeResults, DMatrix_GraphNode, 
                                                    gtfsjs.mSJs[chrchr], gtfexon.GTF_transcript[chrchr], 
                                                    splicechainclass.ClusterCoverage, groupanno.Group_GeneSet, 
                                                    gtfexon.GTF_gene[chrchr], 
                                                    gtfexon.GTF_transcript_strand[chrchr], 
                                                    splicechainclass.HighStrand, chrchr, groupnumber, 
                                                    splicechainclass.HighConClusters,
                                                    FileNo, groupinform.GroupReadFiles,
                                                    SpliceChainSolvent,
                                                    singleexonanno,
                                                    gtfexon.GTF_gene2transcript[chrchr]); 
            // 对每一个文件量化;
            for (int k = 0; k < FileNo; k++) {
                FileFalseNode This_File_False_Node = get_File_False_node(NodeResults.FalseNodeSet,
                                                                        DMatrix_GraphNode.Index2hashname,
                                                                        SpliceChainSolvent.File_HighConClusters,
                                                                        k);
                // 量化;
                // 量化第一步, 准备好所有变量, 初始化;
                IndicateFire InitFirefly = Quantification_initialization_MultiFiles(SpliceChainSolvent.File_ISM[std::to_string(k)], 
                                                            Finally_Annotations, NodeResults.TrueNodeSet,
                                                            This_File_False_Node.ThisFlaseNode, 
                                                            This_File_False_Node.FasleNode2Count,
                                                            DMatrix_GraphNode, groupnumber, 
                                                            SpliceChainSolvent.File_HighConClusters[std::to_string(k)], 
                                                            groupinform.GroupReadSjs, 
                                                            groupinform.GroupReadFiles,
                                                            traceFilePath); //tracemutex
                //量化第二步, EM算法;
                EM_Alg_MultiFiles(Finally_Annotations, InitFirefly, k); // 为了得到所有样本的结果;
                // std::cout << "(=^_^=) Group " << groupinform.GroupIndex << " in file " << k << " completed quantification! (=^_^=) " << std::endl;                     
            }
            AllFile_GeneCounts = get_allfile_genecounts(Finally_Annotations, FileNo, isoformFilePath);
            return AllFile_GeneCounts; // map<gene, map<file,counts>>
        } else {
            // 写到更新的注释中;
            OutputInformation Finally_Annotations = Write_Detection_Transcript2gtf(gtfFilePath, traceFilePath, 
                                                    NodeResults, DMatrix_GraphNode, 
                                                    gtfsjs.mSJs[chrchr], gtfexon.GTF_transcript[chrchr], 
                                                    splicechainclass.ClusterCoverage, groupanno.Group_GeneSet, 
                                                    gtfexon.GTF_gene[chrchr], 
                                                    gtfexon.GTF_transcript_strand[chrchr], 
                                                    splicechainclass.HighStrand, chrchr, groupnumber, 
                                                    splicechainclass.HighConClusters,
                                                    groupinform.GroupReadFiles,
                                                    singleexonanno,
                                                    gtfexon.GTF_gene2transcript[chrchr]);
            // 量化;
            // 量化第一步, 准备好所有变量, 初始化;
            IndicateFire InitFirefly = Quantification_initialization(splicechainclass.ISM, 
                                                        Finally_Annotations, NodeResults.TrueNodeSet,
                                                        NodeResults.FalseNodeSet, 
                                                        DMatrix_GraphNode, groupnumber, 
                                                        splicechainclass.HighConClusters, 
                                                        groupinform.GroupReadSjs,
                                                        groupinform.GroupReadFiles, 
                                                        traceFilePath,
                                                        groupinform.chrName);
            //量化第二步, EM算法;
            AllFile_GeneCounts["0"] = EM_Alg(Finally_Annotations, InitFirefly, isoformFilePath);
            // std::cout << "(=^_^=) Group " << groupinform.GroupIndex << " reads numbers " << groupinform.GroupReadSjs.size()  << " completed quantification! (=^_^=) " << std::endl;  
            return AllFile_GeneCounts; //map<file, <gene,counts>>
        }
    } else {
        return AllFile_GeneCounts;
    }
}


void write_genecount_file(int& FileNo, std::unordered_map<std::string, std::map<std::string, double>>& Single,
                        std::map<std::string, std::map<std::string, double>>& Multi,
                        std::ofstream& geneFilePath){

    if (FileNo == 1) {
        // 只有一个sam文件的情况;
        std::map<std::string, double> thisFileMulti = Multi["0"];
        // 两个文件都有的情况;
        if (Single.size() != 0 && thisFileMulti.size() != 0) {
            // Single多, Multi用来一个个查找;
            if (Single.size() > thisFileMulti.size()) {
                for (const auto & eachG:thisFileMulti) {
                    auto it = Single.find(eachG.first);
                    if (it != Single.end()) {
                        Single[eachG.first]["0"] = Single[eachG.first]["0"] + eachG.second;
                    } else {
                        Single[eachG.first]["0"] = eachG.second;
                    }
                }
                // 开始写入;
                if (geneFilePath.is_open()) {
                    for (const auto& eachG:Single) {
                        std::map<std::string, double> file2count = eachG.second;
                        if (file2count["0"] > 0){
                            std::unique_lock<std::mutex> lock(geneCountMutex);
                            geneFilePath << eachG.first << '\t' << file2count["0"] << '\n';
                        }
                    }
                }                 
            } else {
                // Multi多, Single用来一个个查找;
                for (const auto & eachG:Single) {
                    std::map<std::string, double> file2count = eachG.second;
                    auto it = thisFileMulti.find(eachG.first); 
                    if (it != thisFileMulti.end()) {
                        thisFileMulti[eachG.first] = thisFileMulti[eachG.first] + file2count["0"];
                    } else {
                        thisFileMulti[eachG.first] = file2count["0"];
                    }
                    // 开始写入;
                    if (geneFilePath.is_open()){
                        for (const auto& eachG:thisFileMulti){
                            if (eachG.second > 0){
                                std::unique_lock<std::mutex> lock(geneCountMutex);
                                geneFilePath << eachG.first << '\t' << eachG.second << '\n';
                            }
                        }
                    }
                }
            }
            /* 以上是都有的文件的情况 */
        } else if (thisFileMulti.size() != 0) {
            // 开始写入;
            if (geneFilePath.is_open()){
                for (const auto& eachG:thisFileMulti){
                    if (eachG.second > 0){
                        std::unique_lock<std::mutex> lock(geneCountMutex);
                        geneFilePath << eachG.first << '\t' << eachG.second << '\n';
                    }
                }
            }
            /* 以上是只有多个exon的基因的文件的情况 */
        } else if (Single.size() != 0) {
            // 开始写入;
            if (geneFilePath.is_open()) {
                for (const auto& eachG:Single) {
                    std::map<std::string, double> file2count = eachG.second;
                    if (file2count["0"] > 0) {
                        std::unique_lock<std::mutex> lock(geneCountMutex);
                        geneFilePath << eachG.first << '\t' << file2count["0"] << '\n';
                    }
                }
            }
            /* 以上是只有单个exon的基因的文件的情况 */
        }
    } else {
        if (Single.size() != 0 && Multi.size() != 0) {
            // Single多, Multi用来一个个查找;
            if (Single.size() > Multi.size()) {
                for (const auto& eachG:Multi) {
                    auto it = Single.find(eachG.first);
                    if (it != Single.end()) {
                        for (int i = 0; i < FileNo; i++) {
                            std::string filenumber = std::to_string(i);
                            Single[eachG.first][filenumber] = Single[eachG.first][filenumber] + Multi[eachG.first][filenumber];
                        }
                    } else {
                        for (int i = 0; i < FileNo; i++) {
                            std::string filenumber = std::to_string(i);
                            Single[eachG.first][filenumber] = Multi[eachG.first][filenumber];
                        }                        
                    } 
                }
                // 开始写入文件;
                if (geneFilePath.is_open()) {
                    for (const auto& eachG:Single) {
                        std::map<std::string, double> file2count = eachG.second;
                        std::unique_lock<std::mutex> lock(geneCountMutex);
                        geneFilePath << eachG.first;
                        for (const auto& eachFile:file2count) {
                            geneFilePath << '\t' << eachFile.second;
                        }
                        geneFilePath << '\n';
                    }
                }
            } else {
                for (const auto& eachG:Single) {
                    auto it = Multi.find(eachG.first);
                    if (it != Multi.end()) {
                        for (int i = 0; i < FileNo; i++) {
                            std::string filenumber = std::to_string(i);
                            Multi[eachG.first][filenumber] = Multi[eachG.first][filenumber] + Single[eachG.first][filenumber];
                        }                        
                    } else {
                        for (int i = 0; i < FileNo; i++) {
                            std::string filenumber = std::to_string(i);
                            Multi[eachG.first][filenumber] = Single[eachG.first][filenumber];
                        }    
                    }
                }
                // 开始写入文件;
                if (geneFilePath.is_open()) {
                    for (const auto& eachG:Multi) {
                        std::map<std::string, double> file2count = eachG.second;
                        std::unique_lock<std::mutex> lock(geneCountMutex);
                        geneFilePath << eachG.first;
                        for (const auto& eachFile:file2count) {
                            geneFilePath << '\t' << eachFile.second;
                        }
                        geneFilePath << '\n';
                    }
                }               
            } /* 以上是都有的文件的情况 */
        } else if (Single.size() != 0) {
            // 开始写入文件;
            if (geneFilePath.is_open()) {
                for (const auto& eachG:Single) {
                    std::map<std::string, double> file2count = eachG.second;
                    std::unique_lock<std::mutex> lock(geneCountMutex);
                    geneFilePath << eachG.first;
                    for (const auto& eachFile:file2count) {
                        geneFilePath << '\t' << eachFile.second;
                    }
                    geneFilePath << '\n';
                }
            }            
        } else if (Multi.size() != 0) {
            if (geneFilePath.is_open()) {
                for (const auto& eachG:Multi) {
                    std::map<std::string, double> file2count = eachG.second;
                    std::unique_lock<std::mutex> lock(geneCountMutex);
                    geneFilePath << eachG.first;
                    for (const auto& eachFile:file2count) {
                        geneFilePath << '\t' << eachFile.second;
                    }
                    geneFilePath << '\n';
                }
            }     
        } 
    }
}


// 一个简单的处理函数，可以替换为你的处理逻辑
void processGroup(std::streampos& start, std::streampos& end, 
                  const std::string& sam_file_path, 
                  const int& Sj_supportReadNumber, GTF& gtf_full, 
                  GTFsj& gtf_splice, const int& Mode, 
                  const int& GraphDis,
                  std::ofstream& updatedgtffile, std::ofstream& isoformcountfile,
                  std::ofstream& genecountfile, std::ofstream& tracefile,
                  int& fileno, int& singleEdge) {

    GroupInformation group_information = knowGroupInformation(start, end, sam_file_path, Sj_supportReadNumber);
    std::string chrchr = group_information.chrName;
    std::array<int,2> groupcoverage = group_information.GroupCoverage;
    std::unordered_map<std::string, std::array<int,2>> single_exon_group_annotation;

    std::unordered_map<std::string, std::map<std::string, double>> file_singleexon_gene_number; 
    std::map<std::string, std::map<std::string, double>> file_multiexon_gene_number; 


    if (group_information.GroupSingleExon.size() > 0) {

        single_exon_group_annotation = get_group_single_exon_annotation(gtf_splice.SE[chrchr], groupcoverage);

        std::unordered_map<std::string, std::vector<std::string>> single_exon_with_reads = get_group_single_exon_reads(single_exon_group_annotation, group_information.GroupSingleExon, singleEdge);

        file_singleexon_gene_number = write_single_exon_gtf_trace(fileno, group_information.GroupSingleExon, 
                                    group_information.GroupReadFiles,
                                    single_exon_with_reads, single_exon_group_annotation, 
                                    updatedgtffile, tracefile, isoformcountfile,
                                    gtf_full.GTF_transcript_strand[chrchr], chrchr);
    }
    std::cout << group_information.GroupIndex << " is end! Single reads size is " << group_information.GroupSingleExon.size() << "!"<< std::endl;
    group_information.GroupSingleExon.clear();

    if (group_information.GroupReadSjs.size() > 0) {

        GroupAnnotation group_annotation = get_group_annotation(gtf_splice.mSJsBE[chrchr], gtf_splice.mSJs[chrchr], groupcoverage);

        SpliceChainClass spliceclass = generate_splice_chain_class(Mode, group_information.GroupReadSjs, 
                                                group_information.GroupReadCoverage, 
                                                group_annotation.Group_Annotations, 
                                                group_information.GroupSigns, 
                                                gtf_splice.mSJsBE[chrchr], 
                                                group_information.GroupReadFiles, 
                                                tracefile, Sj_supportReadNumber, chrchr); 
        
        file_multiexon_gene_number = DetectQuant(group_annotation, spliceclass, group_information, GraphDis, gtf_full, gtf_splice, updatedgtffile, isoformcountfile, tracefile, fileno, single_exon_group_annotation); 
    }
    write_genecount_file(fileno, file_singleexon_gene_number, file_multiexon_gene_number, genecountfile);
    std::cout << group_information.GroupIndex << " is end! reads size is " << group_information.GroupReadSjs.size() << "!"<< std::endl;
}


// 主函数;
int main(int argc, char* argv[])
{   
    // 命名空间; 
    int option_index = 0;
    int c;

    std::string samfile_name;
    std::string fastafile_name;
    std::string gtffile_name;
    std::string output_file_name;

    int SJDistance = 18;
    int SJ_support_read_number = 2;
    int Graph_distance = 60;
    int Thread = 10;
    int mode = 0;
    int single_exon_edge = 60;

    while ((c = getopt_long(argc, argv, "s:f:g:o:m:j:n:e:d:t:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 's':
                samfile_name = optarg;
                break;
            case 'f':
                fastafile_name = optarg;
                break;
            case 'g':
                gtffile_name = optarg;
                break;
            case 'o':
                output_file_name = optarg;
                break;
            case 'm':
                mode = std::stoi(optarg);
                break;
            case 'j':
                SJDistance = std::stoi(optarg);  
                break;
            case 'n':
                SJ_support_read_number = std::stoi(optarg);   
                break;
            case 'e':
                single_exon_edge = std::stoi(optarg);   
                break;
            case 'd':
                Graph_distance = std::stoi(optarg);   
                break;
            case 't':
                Thread = std::stoi(optarg);   
                break;
            case 'h':
                print_usage(argv[0]);
                exit(EXIT_FAILURE);              
            case '?':
                std::cerr << "Invalid option" << std::endl;
                exit(EXIT_FAILURE);
            default:
                break;
        }
    }

    if (samfile_name.empty() || fastafile_name.empty() || output_file_name.empty()) {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> sam_file_vec = traverse_sam_file(samfile_name, output_file_name);
    std::vector<std::string> outputFileVec = check_catalog_exist(output_file_name, sam_file_vec);
    std::ofstream gtf_file(outputFileVec[0], std::ios::app);
    std::ofstream isoform_file(outputFileVec[1], std::ios::app);
    std::ofstream gene_file(outputFileVec[2], std::ios::app);
    std::ofstream trace_file(outputFileVec[3], std::ios::app);

    std::cout << "*****" << std::endl;
    std::cout << "SAM file: " << samfile_name << std::endl;
    std::cout << "FASTA file: " << fastafile_name << std::endl;
    std::cout << "GTF file: " << gtffile_name << std::endl;
    std::cout << "Output file: " << output_file_name << std::endl;
    std::cout << "Mode: " << mode << std::endl;
    std::cout << "Single exon boundary: " << single_exon_edge << std::endl;
    std::cout << "SJDistance: " << SJDistance << std::endl;
    std::cout << "SJ_support_read_number: " << SJ_support_read_number << std::endl;
    std::cout << "single_exon_edge : " << single_exon_edge << std::endl;
    std::cout << "Graph_distance: " << Graph_distance << std::endl;
    std::cout << "Thread: " << Thread << std::endl;
    std::cout << "*****" << std::endl;

    std::unordered_map<std::string, std::string> Fasta = Read_fasta_file(fastafile_name);

    FileSplit BroCOLIfile = thread_all_read_sam_files(samfile_name, sam_file_vec, Thread, output_file_name, SJDistance, Fasta);
    Fasta.clear();

    GTF GTF_full = get_gtf_annotation(gtffile_name);

    GTFsj GTF_Splice = get_SJs_SE(GTF_full.GTF_transcript);
    
    std::vector<std::size_t> Group_idx = sort_indexes_e(BroCOLIfile.group_reads_number);
  
    
    ThreadPool BroCOLIpool(Thread);
    
    std::vector<std::future<void>> futures;
    for (const auto& i:Group_idx) {
        futures.emplace_back(BroCOLIpool.enqueue([&, i]() { 
            processGroup(
                BroCOLIfile.reads_pointer[2*i], 
                BroCOLIfile.reads_pointer[2*i+1], 
                BroCOLIfile.readtxt_path, 
                SJ_support_read_number, 
                GTF_full, 
                GTF_Splice,
                mode,
                Graph_distance,
                gtf_file,
                isoform_file,
                gene_file,
                trace_file,
                BroCOLIfile.FileNo,
                single_exon_edge);
            }));
    }

    for (auto& future : futures) {
        future.get();  
    }
    
    gtf_file.close();
    isoform_file.close();
    gene_file.close();
    trace_file.close();
    return 0;
}





