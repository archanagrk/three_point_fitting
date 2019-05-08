
#include "key_struct.h"

/*! \file
 * \struct for three point functions
 */



namespace key_struct
{

  namespace
  {
    //! Error output
    std::ostream& operator<<(std::ostream& os, const Array<int>& d)
    {
      if (d.size() > 0)
      {
	      os << d[0];

	      for(int i=1; i < d.size(); ++i)
	        os << " " << d[i];
      }

      return os;
    }

  }


  //----------------------------------------------------------------------------

  //! Used for output
  std::ostream& operator<<(std::ostream& os, const KeyHadronSUNNPartNPtIrrep_t& op)
  {

    for(int i=0; i < op.irrep.size(); ++i)
    {
      os << " irrep[" << i << "]:  " << op.irrep[i];
    }   

    for(int i=0; i < op.row.size(); ++i)
    {
      os << " row[" << i << "]:  " << op.row[i];
    }   

    for(int i=0; i < op.mom.size(); ++i)
    {
      for(int j=0; j < op.mom[i].size(); ++j )
      {
        os << " mom[" << i << "] [" << j << "]:  " << op.mom[i][j];
      }
    }

    return os;
  }


  //----------------------------------------------------------------------------

  // Read a key
  void read(XMLReader& xml, const std::string& path, KeyHadronSUNNPartNPtIrrep_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "irrep", param.irrep);
    read(paramtop, "row", param.row);
    read(paramtop, "mom", param.mom);

  }


  // Write a key
  void write(XMLWriter& xml, const std::string& path, const KeyHadronSUNNPartNPtIrrep_t& param)
  {
    push(xml, path);

    write(xml, "irrep", param.irrep);
    write(xml, "row", param.row);
    write(xml, "mom", param.mom);


    pop(xml);
  }


  //----------------------------------------------------------------------------


  // Read a key
  void read(BinaryReader& bin, KeyHadronSUNNPartNPtIrrep_t& param)
  {
    for(int i = 0; i < param.irrep.size(); i++){
      readDesc(bin, param.irrep[i]);
    }
    read(bin, param.row);
    read(bin, param.mom);

  }


  // Write a key
  void write(BinaryWriter& bin, const KeyHadronSUNNPartNPtIrrep_t& param)
  {
    for(int i = 0; i < param.irrep.size(); i++){
      writeDesc(bin, param.irrep[i]);
    }
    write(bin, param.row);
    write(bin, param.mom);

  }

} // namespace key_struct






