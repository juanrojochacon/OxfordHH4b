#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath> 

using namespace std;
//copy desired data files into ".cc" files, under names r1.cc, i1.cc, b1.cc for the source, for
//resolved, intermediate and boosted repsectives, with the backgrounds using the following numbers
//e.g. r2.cc is the first resolved background.
//remember to delete all text from file so there are just 5 columns of data!!!
//diHiggs is the signal
//b,c,d,e are the 4 background sources
void signal_significance(string s1="1.cc", string b1="", string b2="", string b3="", string b4="",
			   string s2="", string b5="", string b6="", string b7="", string b8="",
			   string s3="", string b9="", string b10="", string b11="", string b12=""){

  fstream myfile1 (s1), myfile2 (b1), myfile3 (b2), myfile4 (b3), myfile5 (b4);
  fstream myfile6 (s2), myfile7 (b5), myfile8 (b6), myfile9 (b7), myfile10 (b8);
  fstream myfile11 (s3), myfile12 (b9), myfile13 (b10), myfile14 (b11), myfile15 (b12);
  vector<double> vec1 (150), vec2 (150), vec3 (150), vec4(150), vec5(150);
  vector<double> vec6 (150), vec7 (150), vec8 (150), vec9 (150), vec10 (150);
  vector<double> vec11 (150), vec12 (150), vec13 (150), vec14 (150), vec15 (150);
  vector<double> err1 , err2, err3, err4, err5, err6 , err7, err8, err9, err10, err11 , err12, err13, err14, err15;
  vector<double> val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14, val15;  
    
  //Checking the files opened properly 
  if (myfile1.is_open() && myfile2.is_open() && myfile3.is_open() && myfile4.is_open() && myfile5.is_open()
      && myfile6.is_open() && myfile7.is_open() && myfile8.is_open() && myfile9.is_open() && myfile10.is_open()
      && myfile11.is_open() && myfile12.is_open() && myfile13.is_open() && myfile14.is_open(), myfile15.is_open() ){
    cout << "Files opened properly" << endl;}
 

 
    for (int i = 0; i < 150; i++){
      //reading files into vectors
      myfile1 >> vec1[i];
      myfile2 >> vec2[i];
      myfile3 >> vec3[i];
      myfile4 >> vec4[i];
      myfile5 >> vec5[i];
      myfile6 >> vec6[i];
      myfile7 >> vec7[i];
      myfile8 >> vec8[i];
      myfile9 >> vec9[i];
      myfile10 >> vec10[i];
      myfile11 >> vec11[i];
      myfile12 >> vec12[i];
      myfile13 >> vec13[i];
      myfile14 >> vec14[i];
      myfile15 >> vec15[i];
       
      // getting errors from last column of files
      if ((i+1)%5 == 0){
	err1.push_back(vec1[i]);
	err2.push_back(vec2[i]);
	err3.push_back(vec3[i]);
	err4.push_back(vec4[i]);
	err5.push_back(vec5[i]);
	err6.push_back(vec6[i]);
	err7.push_back(vec7[i]);
	err8.push_back(vec8[i]);
	err9.push_back(vec9[i]);
	err10.push_back(vec10[i]);
	err11.push_back(vec11[i]);
	err12.push_back(vec12[i]);
	err13.push_back(vec13[i]);
	err14.push_back(vec14[i]);
	err15.push_back(vec15[i]);
      }
      
     
      // getting values from third column of files
       if((i+3)%5 == 0){
	 val1.push_back(vec1[i]);
	 val2.push_back(vec2[i]);
	 val3.push_back(vec3[i]);
	 val4.push_back(vec4[i]);
	 val5.push_back(vec5[i]);
	 val6.push_back(vec6[i]);
	 val7.push_back(vec7[i]);
	 val8.push_back(vec8[i]);
	 val9.push_back(vec9[i]);
	 val10.push_back(vec10[i]);
	 val11.push_back(vec11[i]);
	 val12.push_back(vec12[i]);
	 val13.push_back(vec13[i]);
	 val14.push_back(vec14[i]);
	 val15.push_back(vec15[i]);

       }
    }

   
      
       //calculating total signal and total  backgrounds
       double sig = 0, back1 = 0, back2 = 0, back3 = 0, back4 = 0;
      
       for (unsigned int l = 0; l < val1.size(); l++){
	 sig += val1[l] + val6[l] + val11[l];
	 back1 += val2[l] + val7[l] + val12[l];
	 back2 += val3[l] + val8[l] + val13[l];
	 back3 += val4[l] + val9[l] + val14[l];
	 back4 += val5[l] + val10[l] + val15[l];
       }
       cout << "Total Signal (Sum of Weights for diHiggs) = " << sig << endl;
	
       
  // Cross sections obtained via xsection.cc, and pheno paper for xsec_sig, measured in pb^-1
       double xsec_sig = 6.2 * pow(10,-3), xsec_back1 = 1.1 * pow(10,3),
	 xsec_back2 = 2.7 * pow(10,5) , xsec_back3 = 9.9 * pow(10,6), xsec_back4 = 5.7 * pow(10,2);
      
       // number of trials before cuts obtained via xsection.cc
       // HAVE JUST LEFT "ns = 1" BUT THIS NEEDS TO BE CHANGED
       double ns = 1, nb1 = 6.14405*pow(10,9), nb2 = 7.41471*pow(10,9), nb3 = 9.27505*pow(10,9) , nb4 = 4.95184*pow(10,7);

       // calculating effective cross sections of signal and each background in signal region in pb^-1
       double x_s = sig/ns;
       double x_b1 = back1/nb1, x_b2 = back2/nb2, x_b3 = back3/nb3, x_b4 = back4/nb4;
   
       // calculating effective event numbers of signal and each background in signal region
       double L_ref = 3000000; //pb^-1
       double n_s = L_ref*x_s, n_b1 = L_ref*x_b1, n_b2 = L_ref*x_b2, n_b3 = L_ref*x_b3, n_b4 = L_ref*x_b4;    
       
       cout << "Sum of trials before cuts for 2b2j = " << nb1 << endl;
       cout << "Sum of trials before cuts for 4b = " << nb2 << endl;
       cout << "Sum of trials before cuts for 4j = " << nb3 << endl;
       cout << "Sum of trials before cuts for ttbar->X = " << nb4 << endl;

       cout << "Sum of weights after cuts for 2b2j = " << back1 << endl;
       cout << "Sum of weights after cuts for 4b = " << back2 << endl;
       cout << "Sum of weights after cuts for 4j = " << back3 << endl;
       cout << "Sum of weights after cuts for ttbar->X = " << back4 << endl;

      
       cout << "Effective number of 2b2j events = " << n_b1 << endl;
       cout << "Effective number of 4b events = " << n_b2 << endl;
       cout << "Effective number of 4j events = " << n_b3 << endl;
       cout << "Effective number of ttbar->X events = " << n_b4 << endl;


    myfile1.close();
    myfile2.close();
    myfile3.close();
    myfile4.close();
    myfile5.close();
    
    
    
    
}

int main(){

  signal_significance("1.cc","2.cc","3.cc","4.cc","5.cc","6.cc","7.cc","8.cc","9.cc","10.cc","11.cc","12.cc","13.cc","14.cc","15.cc");

  
 return 0;
}
