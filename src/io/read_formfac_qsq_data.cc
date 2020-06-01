#include "read_formfac_qsq_data.h"

map<double , ENSEM::EnsemReal> readlist::make_ff_qsq(const string& listfile){

    map<double , ENSEM::EnsemReal> out;

    string line, buf; ifstream list(listfile.c_str());
    if(!(list.is_open())){ std::cerr << __func__ << ": error reading " << listfile << std::endl; exit(1); }
    
    while( !(list.eof()) ){
      if( !getline(list, line) ) break;
      bool active = true;
      if(line.find("(*)") != string::npos){ active = false; }

      /* check for deactivate flag */
      if(active){
    
        vector<string> tokens; stringstream ss(line); while(ss >> buf){tokens.push_back(buf);} /* split by whitespace */
        
        if(tokens.size() < 2){ cerr << "not enough data in line : " << line << endl; exit(1);}


        double qsq; string ff_file; ENSEM::EnsemReal ff;
      
        { qsq = stod(tokens[0].c_str()); }

        ff_file = tokens[1].c_str();
        read(ff_file,ff);

        if( out.find(qsq) != out.end() ){ cerr << "duplicate qsq in input file :: " << qsq << ", exiting" << endl; exit(1); }  

        out.insert(make_pair(qsq,ff));

    }
  }
  list.close();
  
  return out;
};
