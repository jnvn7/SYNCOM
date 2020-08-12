// SYNCOM - A Nonlinear Synthetic Rope Numerical Computation Software
//
// Copyright (c) 2020 Jessica Nguyen <nvnguyen@umass.edu>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
////////////////////////////////////////////////////////////////////////////////


#ifndef readIn_api_h
#define readIn_api_h

#include "setting.h"
#include "error.h"
#include "rapidxml-1.13/rapidxml.hpp"
#include "rapidxml-1.13/rapidxml_print.hpp"
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include <regex>

namespace rope {
    // Group of functions used to read the input file.

    class ReadIn
    {

    public:

        // Read input data from main input file.
        ReadIn(void) {};

        // Read input data from the users' defined XML file. 
        // Data includes material properties of ropes and applied stress (or strain);
        ErrorCode readIn_data(Setting& setting);

    private:
        
        // Used when reading main input data file.
        int readInput(std::vector<std::vector<double>>& dataIn, const std::string data_file_name, const int header_rows);
        int check_file_existence(const std::string file_name);
        int check_availability(const rapidxml::xml_node<>* node,
            const std::vector<std::string>& names);
        int check_is_number(const rapidxml::xml_node<>* node,
            const std::vector<std::string>& names);
        int extract_multi_number(const std::string token,
            std::vector<std::string>& number_string);
        bool is_number(const std::string& token);
        bool is_integer(const std::string& token);
    };

} // End of namespace rope.

#endif // readIn_api_h
