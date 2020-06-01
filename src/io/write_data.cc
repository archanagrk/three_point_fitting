

#include "write_data.h"

//************************************************************************
void print_output(string& path, string& name, prefactor& pf, fit_three_point_output& output_rl, fit_three_point_output& output_im,
        string& Zfile, map<string, vector< pair< pair<ENSEM::EnsemReal, ENSEM::EnsemReal> , ENSEM::EnsemReal >>>& fqsq,
        map<string, vector<ENSEM::EnsemReal> >& favg, map<string, vector<ENSEM::EnsemReal> >& fmean)

{


    if(output_rl.success){


      // /* write the Q^2 vs F(Q^2) output */
      {

        ENSEM::EnsemReal fq2 = output_rl.F; 

        {
          ostringstream outfile; outfile << path << "Real/";
          SEMBLEIO::makeDirectoryPath(outfile.str());

          
          outfile << name << "_FReal.jack"; 
          write(outfile.str(), fq2);
        }

      }

      /* write log to file */
      {
        stringstream s; s << path << "Real/" << name << "_three_pt_fit_real.log"; 
        ofstream out; out.open(s.str().c_str());
        out << output_rl.fit_summary;
        out.close();
      }
      
      /* write mass ensem file */
      {  ostringstream outfile; outfile << path << "Real/" <<  name << "_three_pt_fit_real.jack";  write(outfile.str(), output_rl.F );  }

      /* write mass variations to file */
      {
        stringstream s; s << path << "Real/" <<  name << "_three_pt_fit_real.syst"; 
        ofstream out; out.open(s.str().c_str());
        out << output_rl.F_fit_variation;
        out.close();
      }

      // /* write plot data to file */
      {
        stringstream s; s << path << "Real/" <<  name << "_three_pt_fit_real.plot"; 
        ofstream out; out.open(s.str().c_str());
        out << "## name= " << name << endl;
        out << output_rl.plot_data;
        out.close();
      }

    }


    if(output_im.success){

      // /* write the Q^2 vs F(Q^2) output */
      {

        ENSEM::EnsemReal fq2 = output_im.F;//*Zv*Zs*CurrFac;

        {
          ostringstream outfile; outfile << path << "Imag/";
          SEMBLEIO::makeDirectoryPath(outfile.str());

          outfile << name << "_FImag.jack"; 
          write(outfile.str(), fq2);
        }
        
      }

      /* write log to file */
      {
        stringstream s; s << path << "Imag/" << name << "_three_pt_fit_imag.log"; 
        ofstream out; out.open(s.str().c_str());
        out << output_im.fit_summary;
        out.close();
      }
      
      /* write mass ensem file */
      {  ostringstream outfile; outfile << path << "Imag/" <<  name << "_three_pt_fit_imag.jack";  write(outfile.str(), output_im.F );  }

      /* write mass variations to file */
      {
        stringstream s; s << path << "Imag/" << name << "_three_pt_fit_imag.syst"; 
        ofstream out; out.open(s.str().c_str());
        out << output_im.F_fit_variation;
        out.close();
      }

      // /* write plot data to file */
      {
        stringstream s; s << path << "Imag/" << name << "_three_pt_fit_imag.plot"; 
        ofstream out; out.open(s.str().c_str());
        out << "## name= " << name << endl;
        out << output_im.plot_data;
        out.close();
      }

    }

    //=============================================================================
    //======= DECIDE IF THE IMAGINARY PART IS SIGNIFICANT AND OUTPUT STUFF ======== 

    ENSEM::EnsemReal fq2;


    bool out_rl = false;
    bool out_im = false;
    bool write_out = true;
    pair<double,double>  const_rl, const_im;
    
    if(output_rl.success){
      const_rl = mean_err(output_rl.F);
      if( (abs(const_rl.first) - abs(3.*(const_rl.second)) ) > 0. ){out_rl = true;}
    }
    if(output_im.success){
      const_im = mean_err(output_im.F);
      if( (abs(const_im.first) - abs(3.*(const_im.second)) ) > 0. ){out_im = true;}
    }
    

    if(output_rl.success && output_im.success){

      EnsemReal ratio_ensem = rescaleEnsemUp(  atan( rescaleEnsemDown(output_im.F)/rescaleEnsemDown(output_rl.F) ));
      pair<double,double> ratio = mean_err(ratio_ensem);     

      if( (0.25 > (abs(ratio.first) + ratio.second)) && ( (abs(ratio.first) - ratio.second) > -0.25) ){


        ostringstream outfile; outfile << path << "fit.log";
        ofstream file; // out file stream
        file.open(outfile.str());
        file << "The Arg(F) is: " << ratio.first << " ( "  << ratio.second << " ) " << endl;
        file << "the real part is: " << const_rl.first << " ( " << const_rl.second << " ) " << endl;
        file << "the imag part is: " << const_im.first << " ( " << const_im.second << " ) " << endl;
        file << "****************************************************************" << endl << endl;
        if(out_rl){
          file << "The value is real with the value of F: " << const_rl.first << " ( " << const_rl.second << " ) " << endl;

          file.close();
              
          fq2 = output_rl.F;

          {
            ostringstream outfile; outfile << path << name << "_F.jack"; 
            write(outfile.str(), fq2);
          }

          if(Zfile == "true"?true:false)
          {

            ENSEM::EnsemReal Zv = toScalar(1.0)/output_rl.F; 

            {
              ostringstream outfile; outfile << path << name << "_Zv.jack"; 
              write(outfile.str(), Zv);

            }        

          }          
          
        }

        else{file << "Decided the value is statistically consistent with zero" << endl; file.close();  write_out = false;  }

        
      }    


      else{

        ostringstream outfile; outfile << path << "fit.log";
        ofstream file; // out file stream
        file.open(outfile.str());
        file << "The Arg(F) is: " << ratio.first << " ( "  << ratio.second << " ) " << endl;
        file << "the real part is: " << const_rl.first << " ( " << const_rl.second << " ) " << endl;
        file << "the imag part is: " << const_im.first << " ( " << const_im.second << " ) " << endl;
        file << "****************************************************************" << endl << endl;

        if(out_im){
          file << "The value is imag with the value of F: " << const_im.first << " ( " << const_im.second << " ) " << endl;

          file.close();
              
          fq2 = output_im.F;

          {
            ostringstream outfile; outfile << path << name << "_F.jack"; 
            write(outfile.str(), fq2);
          }

          if(Zfile == "true"?true:false)
          {

            ENSEM::EnsemReal Zv = toScalar(1.0)/output_im.F; 

            {
              ostringstream outfile; outfile << path << name << "_Zv.jack"; 
              write(outfile.str(), Zv);

            }        

          }
        }

        else{file << "Decided the value is statistically consistent with zero" << endl; file.close(); write_out = false; }

      }

    }

    else if(output_rl.success){
      ostringstream outfile; outfile << path << "fit.log";
      ofstream file; // out file stream
      file.open(outfile.str());

      
      file << "The fits to the imaginary part failed" << endl;
      file << "The value of the real part is: " << const_rl.first << " ( " << const_rl.second << " ) " << endl;
      file << "****************************************************************" << endl << endl;
      if(out_rl){
        file << "The value is real with the value of F: " << const_rl.first << " ( " << const_rl.second << " ) " << endl; 

        file.close();
          
        fq2 = output_rl.F;

        {
          ostringstream outfile; outfile << path << name << "_F.jack"; 
          write(outfile.str(), fq2);
        }

        if(Zfile == "true"?true:false)
        {

          ENSEM::EnsemReal Zv = toScalar(1.0)/output_rl.F; 

          {
            ostringstream outfile; outfile << path << name << "_Zv.jack"; 
            write(outfile.str(), Zv);

          }        

        }
      }
      else{file << "Decided the value is statistically consistent with zero" << endl; file.close(); write_out = false; }

      cout << "All fits to the imaginary part failed " << endl;
    }

    else if(output_im.success){
      ostringstream outfile; outfile << path << "fit.log";
      ofstream file; // out file stream
      file.open(outfile.str());

      file << "The fits to the real part failed" << endl;
      file << "The value of the imag part is: " << const_im.first << " ( " << const_im.second << " ) " << endl;
      file << "****************************************************************" << endl << endl;
      if(out_im){
        file << "The value is imag with the value of F: " <<  const_im.first << " ( " << const_im.second << " ) " << endl;
        file.close();
            
        fq2 = output_im.F;

        {
          ostringstream outfile; outfile << path << name << "_F.jack"; 
          write(outfile.str(), fq2);
        }

        if(Zfile == "true"?true:false)
        {

          ENSEM::EnsemReal Zv = toScalar(1.0)/output_im.F; 

          {
            ostringstream outfile; outfile << path << name << "_Zv.jack"; 
            write(outfile.str(), Zv);

          }        

        }
        
        cout << "All fits to the real part failed " << endl;
      }

      else{file << "Decided the value is statistically consistent with zero" << endl; file.close(); write_out = false; }
    }

    else{ cout << "All fits to the real part and the imaginary part failed" << endl; write_out = false; }



    //=============================
    //======= OUTPUT STUFF ======== 

    if(write_out){

      SEMBLEIO::makeDirectoryPath(path + "../jackfiles/"); 

      {
        ostringstream outfile; outfile << path << name << "_Q2.jack"; 
        write(outfile.str(), pf.qsq);
      }

      {
        ostringstream outfile; outfile << path << "../jackfiles/" << "Q2.jack"; 
        write(outfile.str(), pf.qsq);
      }

      {
        ostringstream outfile; outfile << path << name << "_Estar.jack"; 
        write(outfile.str(), pf.ei);
      }

      if(toDouble(mean(fq2)) == toDouble(mean(output_rl.F))){
        if(fqsq.find(name) == fqsq.end() ){

          vector<pair<pair< ENSEM::EnsemReal, ENSEM::EnsemReal>, ENSEM::EnsemReal>> fqsq_val;
          fqsq_val.push_back(make_pair(make_pair(pf.qsq, pf.ei),fq2));
          fqsq.insert(make_pair(name, fqsq_val)); 

          }

        else{ fqsq.find(name)->second.push_back(make_pair(make_pair(pf.qsq, pf.ei),fq2)); }

        stringstream s; s <<  std::fixed << std::setprecision(6) << toDouble(mean(pf.qsq));
        string q2 = s.str();

        stringstream ss; ss <<  std::fixed << std::setprecision(6) << toDouble(mean(pf.ei));
        string ei = ss.str();

        if(favg.find(q2) == favg.end()){ vector<ENSEM::EnsemReal> val; val.push_back(fq2); favg.insert( make_pair(q2,val) );}
        else{favg.find(q2)->second.push_back(fq2);}

        if(fmean.find(ei) == fmean.end()){ vector<ENSEM::EnsemReal> val; val.push_back(fq2); fmean.insert( make_pair(ei,val) );}
        else{fmean.find(ei)->second.push_back(fq2);}
      }
    }

};



//************************************************************************
void print_extra_data(string& xmlini, string& Zfile, map<string, vector< pair< pair<ENSEM::EnsemReal, ENSEM::EnsemReal> , ENSEM::EnsemReal >>>& fqsq,
        map<string, vector<ENSEM::EnsemReal> >& favg)

  {

    std::string path = SEMBLEIO::getPath();
    std::string outdir = path + "output/";
    SEMBLEIO::makeDirectoryPath(outdir);


    {
      stringstream s; s << outdir <<  xmlini << "_FQ2.out";         
      ofstream out; out.open(s.str().c_str());
      out << "## Q  Q_E  F  F_E" << endl << endl;
      for(auto it = fqsq.begin(); it != fqsq.end(); it++ )
      {

        out << "## irrep=" << it->first << endl;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){  
          pair<double, double> q_out = mean_err(it1->first.first);
          pair<double, double> f_out = mean_err(it1->second);
          out << q_out.first << "  " << q_out.second << "  " <<  f_out.first <<  "  " << f_out.second << endl;}
          out << endl << endl;
      }

      out.close();  
    }

    {
      stringstream s; s << outdir <<  xmlini << "_ZvQ2.out";         
      ofstream out; out.open(s.str().c_str());
      out << "## Q  Q_E  F  F_E" << endl << endl;
      for(auto it = fqsq.begin(); it != fqsq.end(); it++ )
      {

        out << "## irrep=" << it->first << endl;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){  
          pair<double, double> q_out = mean_err(it1->first.first);
          pair<double, double> f_out = mean_err(toScalar(1)/it1->second);
          out << q_out.first << "  " << q_out.second << "  " <<  f_out.first <<  "  " << f_out.second << endl;}
          out << endl << endl;
      }

      out.close();  
    }

    {
      stringstream s; s << outdir <<  xmlini << "_FMv.out";         
      ofstream out; out.open(s.str().c_str());
      out << "## E  E_E  F  F_E" << endl << endl;
      for(auto it = fqsq.begin(); it != fqsq.end(); it++ )
      {

        out << "## irrep=" << it->first << endl;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){ 
          pair<double, double> e_out = mean_err(it1->first.second);
          pair<double, double> f_out = mean_err(it1->second);
          out << e_out.first << "  " << e_out.second << "  " <<  f_out.first <<  "  " << f_out.second << endl;}
          out << endl << endl;
      }

      out.close();  
    }

    std::string plotdir = path + "plots/";
    SEMBLEIO::makeDirectoryPath(plotdir); 



    if(Zfile == "true"?true:false){
      stringstream ff; ff << plotdir <<  xmlini << "_FF.plot";   
      stringstream Zvf; Zvf << plotdir <<  xmlini << "_Zv.plot";   
            
      ofstream write_ff; write_ff.open(ff.str().c_str());
      write_ff << "## Q2        F        F_E" << endl << endl;

      ofstream write_zvf; write_zvf.open(Zvf.str().c_str());
      write_zvf << "## Q2        Zv        Zv_E" << endl << endl;

      for(auto it = favg.begin(); it != favg.end(); it++ ){
        ENSEM::EnsemReal avg = rescaleEnsemDown(it->second.at(0));
        for(auto it1 = it->second.begin()+1; it1 != it->second.end(); it1++){avg += rescaleEnsemDown(*it1);}

        avg = rescaleEnsemUp(avg)/toScalar(it->second.size());

        stringstream s; s << path << "Q2_" << it->first << "/jackfiles/FF.jack";
        write(s.str(),avg);
        stringstream ss; ss << path << "Q2_" << it->first << "/jackfiles/Zv.jack"; ENSEM::EnsemReal Zavg = toScalar(1)/avg; write(ss.str(), Zavg);  

        pair<double, double> f_avg =  mean_err(avg); 
        pair<double, double> Z_avg =  mean_err(toScalar(1)/avg);    

        write_ff << it->first << " " << f_avg.first << " " <<  f_avg.second << endl;
        write_zvf << it->first << " " << Z_avg.first << " " <<  Z_avg.second << endl;

      }

      write_ff.close(); 
      write_zvf.close();


      
    }


    else{
      stringstream ff; ff << plotdir <<  xmlini << "_FF.plot";   
  
            
      ofstream write_ff; write_ff.open(ff.str().c_str());
      write_ff << "## Q2        F        F_E" << endl << endl;


      for(auto it = favg.begin(); it != favg.end(); it++ ){
        ENSEM::EnsemReal avg = rescaleEnsemDown(it->second.at(0));

        for(auto it1 = it->second.begin()+1; it1 != it->second.end(); it1++){avg +=  rescaleEnsemDown(*it1);}

        avg = rescaleEnsemUp(avg)/toScalar(it->second.size());
        stringstream s; s << path << "Q2_" << it->first << "/FF.jack";


        write(s.str(),avg);

        pair<double, double> f_avg =  mean_err(avg);     
        write_ff << it->first << " " << f_avg.first << " " <<  f_avg.second << endl;


      }

      write_ff.close(); 
    }


  }
