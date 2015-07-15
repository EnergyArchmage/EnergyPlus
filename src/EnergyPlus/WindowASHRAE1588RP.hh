
// EnergyPlus Headers
#include <EnergyPlus.hh>

// JSON Header
#include <json/json.h>

namespace EnergyPlus {

namespace WindowASHRAE1588RP {

void
CreateASHRAE1588RPConstructions( int & ConstrNum, bool & ErrorsFound );

Json::Value read_1588_database(std::string file_path);

void search_database_keys_for_input(const std::string &construction_name, const std::string &field_name, const std::vector< std::string > &keys, const std::string &input, bool &ErrorsFound);

void create_dummy_variables();

void remove_dummy_variables();

void calc_window_performance(Real64 T_in, Real64 T_out, Real64 v_ws, Real64 I_s);

} // WindowASHRAE1588RP

} // EnergyPlus
