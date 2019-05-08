  //============================
  //===== INTERFACE TO EDB =====

  Register all the handlers
  bool success = edb_reader::registerAll();

  for(int i = 0; i < num; i++ ){ 

    // Read the dbtype from the metadata of the first file
    int SIZE = 3;
    char* edb_in[] = {"sample.edb","keysxml","all"};

    if (SIZE <= 2){return edb_reader::usage(SIZE, edb_in);}

    string dbtype = getDBType(edb_in[0]);

    // Call the appropriate handler
    int ret;
    ret = TheKeyHandlerObjFuncMap::Instance().callFunction(dbtype, SIZE, edb_in);

    char* dat_out[] = {"sample.edb","get","all"};

    // Parse the arguments
    if (SIZE <= 2){return edb_reader::usage(SIZE, dat_out);}

    dbtype = getDBType(dat_out[0]);
    ret = TheKeyHandlerObjFuncMap::Instance().callFunction(dbtype, SIZE, dat_out);

  }

  //============================
  //===== LOAD THE DATA =====

  /* Generates vector of data with the size  = number of Dt. Each data[i] has data for a particular Dt and all the Dts less than it.
   First ensemble fits are done on the smallest Dt. The range is fixed and then moves to the next Dt and does a combined fit with the fixed range
   smaller Dt and the curret Dt to fix the range of the current Dt. Then finally data[num - 1] contains all the data the the output of this fit is used.
   reduces the complexity from n^m to n*m */


  vector<pair<int,int>> x; vector<EnsemReal> y_ensem; EnsemVectorReal tmp;
  vector<pair<int,int>> x_all; vector<EnsemReal> y_ensem_all;
  vector<Data> data(num+1);


  for(int i = 0; i < num; i++ ){

    read(filename[i], tmp);

    for(int t = 0; t < tmp.numElem(); t++){

      x.push_back(make_pair(dt_v[i],t));
      y_ensem.push_back( peekObs(tmp, t) );

      x_all.push_back(make_pair(dt_v[i],t));
      y_ensem_all.push_back( peekObs(tmp, t) );

    }

    make_pair_int_abscissa_data(x, y_ensem, data[i] );
    x.clear(); y_ensem.clear();

  }

  make_pair_int_abscissa_data(x_all, y_ensem_all, data[num] );
  cout << "loaded data::" << endl << data[num].print_data() << endl;
