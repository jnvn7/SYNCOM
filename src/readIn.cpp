// SYNCOM - A Nonlinear Synthetic Rope Numerical Computation Software
//
// Copyright 2020 Jessica Nguyen <nvnguyen@umass.edu>
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


#include "readIn.h"

using namespace std;
using namespace rapidxml;

namespace rope {

    ////////////////////////////////////////////////////////////////////////////////
    /// Initialization of the setting.
    ////////////////////////////////////////////////////////////////////////////////
    ErrorCode ReadIn::readIn_data(Setting& setting)
    {
        // Get the absolute path to Setting.xml file.
        std::string file_name = setting.setting_folder + setting.setting_file;
        
        ////////////////////////////////////////////////////////////////////////////
        // Setting file.
        ////////////////////////////////////////////////////////////////////////////
        int flag = check_file_existence(file_name); 
        if (flag) 
            return ErrorCode::SETTING_FILE_NONEXISTENT;
        
        ifstream file(file_name);
        std::stringstream buffer;
        buffer << file.rdbuf();
        file.close();

        try
        {
            xml_document<> set_doc;
            std::string content(buffer.str());
            set_doc.parse<0>(&content[0]);
        }
        catch (const rapidxml::parse_error& e)
        {
            return ErrorCode::SETTING_FILE_ERROR_PARSE;
        }

        xml_document<> set_doc;
        std::string content(buffer.str());
        set_doc.parse<0>(&content[0]);

        xml_node<>* root_node;
        root_node = set_doc.first_node("settings");
        if (root_node == 0) 
            return ErrorCode::SETTING_FILE_NO_SETTING_NODE; 
       
        ////////////////////////////////////////////////////////////////////////////
        // Check Data (stress or strain) input file.
        ////////////////////////////////////////////////////////////////////////////
        // Get the data input file.
        xml_node<>* child_node = root_node->first_node("data_input_file");
        if (child_node == 0)
            return ErrorCode::SETTING_FILE_NO_INPUT_DATA_FILE;

        // Check whether the file exists.
        std::string relative_path = child_node->value();
        setting.input_data_path = setting.setting_folder + relative_path;

        flag = check_file_existence(setting.input_data_path);
        if (flag) return ErrorCode::INPUT_DATA_FILE_NONEXISTENT;

        size_t work_folder_index = setting.input_data_path.find_last_of("/\\");
        setting.work_folder = setting.input_data_path.substr(0, work_folder_index + 1); 
        setting.output_filename = setting.work_folder + "SYNCOM_Ouput";
        setting.log_filename = setting.work_folder + "SYNCOM_Log.txt";

        ////////////////////////////////////////////////////////////////////////////
        // Check Module Selection and Material Properties Input Data.
        ////////////////////////////////////////////////////////////////////////////
        // Module selection;
        child_node = root_node->first_node("module");
        if (child_node == 0)
            return ErrorCode::SETTING_FILE_NO_MODULE_SELECTION;

        if (!is_integer(child_node->value()))
            return ErrorCode::SETTING_FILE_BAD_MODULE_SELECTION;
        
        setting.module = stoi(child_node->value());
     
        // Material properties;
        child_node = root_node->first_node("material_props");
        if (child_node == 0)
            return ErrorCode::SETTING_FILE_NO_MATERIAL_PROPERTIES;

        std::vector<std::string> names;
        names.push_back("sigma_yield0");
        names.push_back("MBL");
        names.push_back("Do");
        names.push_back("Dn");

        if (!check_availability(child_node, names))
            return ErrorCode::SETTING_FILE_INCOMPLETE_MATERIAL_PROPERTIES;
        else
        {
            names.erase(names.end()-1);
            if (!check_is_number(child_node, names))
                return ErrorCode::SETTING_FILE_NAN_MATERIAL_PROPERTIES;
        }

        names.clear();

        setting.material_props->sigma_yield0 =
            stod(child_node->first_node("sigma_yield0")->value());

        setting.material_props->MBL =
            stod(child_node->first_node("MBL")->value());

        setting.material_props->Do =
            stod(child_node->first_node("Do")->value()); 
        
        /// Read in Dn's values;
        std::vector<std::string> DnString;
        flag = extract_multi_number(child_node->first_node("Dn")->value(), DnString);
       
        if (DnString.size() <= 0)
            return ErrorCode::SETTING_FILE_BAD_DN_VALUES;
        else
        {
            setting.material_props->Dn.resize(DnString.size());
            setting.material_props->lamdaN.resize(DnString.size());
            for (int i = 0; i < DnString.size(); i++) {
                setting.material_props->Dn[i] = stod(DnString[i]);
                if (stod(DnString[i]) < 0)
                    return ErrorCode::BAD_MATERIAL_PROPERTIES_INPUT;

                setting.material_props->lamdaN[i] = pow(10,-i);
                setting.material_props->sumDn += setting.material_props->Dn[i];
            }
        }
        
        /// Read in a0, g0, g1, g2, Ep, np, H's values;
        std::vector<std::string> a0String;
        std::vector<std::string> g0String;
        std::vector<std::string> g1String;
        std::vector<std::string> g2String;
        std::vector<std::string> EpString;
        std::vector<std::string> npString;
        std::vector<std::string> HString;

        extract_multi_number(child_node->first_node("a0")->value(), a0String);
        extract_multi_number(child_node->first_node("g0")->value(), g0String);
        extract_multi_number(child_node->first_node("g1")->value(), g1String);
        extract_multi_number(child_node->first_node("g2")->value(), g2String);
        extract_multi_number(child_node->first_node("Ep")->value(), EpString);
        extract_multi_number(child_node->first_node("np")->value(), npString);
        extract_multi_number(child_node->first_node("H")->value(), HString);
        
        if (a0String.size() <= 0 || g0String.size() <= 0 || g1String.size() <= 0
            || g2String.size() <= 0 || EpString.size() <= 0 || npString.size() <= 0
            || HString.size() <= 0) 
            return ErrorCode::SETTING_FILE_NAN_MATERIAL_PROPERTIES;
        else
        {
            setting.material_props->a0Coefs.resize(a0String.size());
            setting.material_props->g0Coefs.resize(g0String.size());
            setting.material_props->g1Coefs.resize(g1String.size());
            setting.material_props->g2Coefs.resize(g2String.size());
            setting.material_props->EpCoefs.resize(EpString.size());
            setting.material_props->npCoefs.resize(npString.size());
            setting.material_props->H_vpCoefs.resize(HString.size());

            for (size_t i = 0; i < a0String.size(); i++) {
                setting.material_props->a0Coefs[i] = stod(a0String[i]);
            }
            for (size_t i = 0; i < g0String.size(); i++) {
                setting.material_props->g0Coefs[i] = stod(g0String[i]);
            }
            for (size_t i = 0; i < g1String.size(); i++) {
                setting.material_props->g1Coefs[i] = stod(g1String[i]);
            }
            for (size_t i = 0; i < g2String.size(); i++) {
                setting.material_props->g2Coefs[i] = stod(g2String[i]);
            }
            for (size_t i = 0; i < EpString.size(); i++) {
                setting.material_props->EpCoefs[i] = stod(EpString[i]);
            }
            for (size_t i = 0; i < npString.size(); i++) {
                setting.material_props->npCoefs[i] = stod(npString[i]);
            }
            for (size_t i = 0; i < HString.size(); i++) {
                setting.material_props->H_vpCoefs[i] = stod(HString[i]);
            }
        }
        
        ////////////////////////////////////////////////////////////////////////////
        // Check Material Properties Input Data.
        ////////////////////////////////////////////////////////////////////////////
        child_node = root_node->first_node("numerical_setting");
        if (child_node == 0)
            return ErrorCode::SETTING_FILE_NO_MATERIAL_PROPERTIES;

        names.push_back("dt");
        names.push_back("limit");
        names.push_back("tol");

        if (!check_availability(child_node, names))
            return ErrorCode::SETTING_FILE_INCOMPLETE_MATERIAL_PROPERTIES;
        else
        {
            if (!check_is_number(child_node, names))
                return ErrorCode::SETTING_FILE_NAN_MATERIAL_PROPERTIES;
        }

        setting.dt = stod(child_node->first_node("dt")->value());

        setting.limit = stoi(child_node->first_node("limit")->value());

        setting.tol = stod(child_node->first_node("tol")->value()); 
    
        ////////////////////////////////////////////////////////////////////////////
        // Read input stress (strain) data from provided file.
        ////////////////////////////////////////////////////////////////////////////
        flag = readInput(setting.dataIn, setting.input_data_path, 1);

        switch (flag) {
        case 1: 
            return ErrorCode::FAIL_TO_OPEN_INPUT_FILE;
            break;
        case 2:
            case 3: 
                return ErrorCode::WRONG_INPUT_FILE_FORMAT;
                break;
        case 4: 
            return ErrorCode::NAN_INPUT_DATA;
            break;
        }

        return ErrorCode::SUCCESS;
    }
    ////////////////////////////////////////////////////////////////////////////////
    /// Read data matrix with header lines for stress/strain users' input data.
    ////////////////////////////////////////////////////////////////////////////////
    int ReadIn::readInput(std::vector<std::vector<double>>& dataIn, const std::string data_file_name, const int header_rows)
    {
        string line, cell;
        int i_line = 0, i_cell;
        int expected_cols = 1;
        ifstream data_file(data_file_name);

        if (!data_file.good() || (data_file_name.back() == '/'
            || data_file_name.back() == '\\'
            || data_file_name.back() == '.'))
        {
            data_file.close();
            return 1; // Open file failed.
        }
        else
        {
            // Get number of lines.
            int i_line = 0;
            while (std::getline(data_file, line))
                i_line++;

            int n_points = i_line - 1;
            dataIn.resize(n_points, vector<double>(expected_cols));
            
            data_file.close();
            data_file.open(data_file_name);

            i_line = 0;
            while (std::getline(data_file, line))
            {
                std::stringstream line_stream(line);
                if (i_line >= header_rows)
                {
                    i_cell = 0;
                    while (line_stream >> cell)
                    {
                        if (!(i_cell < expected_cols))
                            return 3; // More than expected columns found.
                        else if (!is_number(cell))
                            return 4; // NaN found in data.
                        else
                        {
                            dataIn[i_line - header_rows][i_cell] = stod(cell);
                        }
                        i_cell++;
                    }
                    if (i_cell != expected_cols)
                        return 2;
                }
                i_line++;
            }
            data_file.close();
            return 0; // Success.
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Check whether file exist.
    ////////////////////////////////////////////////////////////////////////////////
    int ReadIn::check_file_existence(const std::string file_name)
    {
        ifstream file(file_name);
        if (!file.good() || (file_name.back() == '/' || file_name.back() == '\\'
            || file_name.back() == '.'))
            return 1; // Open file failed.
        else
        {
            file.close();
            return 0;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Check whether expected subnodes or attributes are available.
    ////////////////////////////////////////////////////////////////////////////////
    int ReadIn::check_availability(const rapidxml::xml_node<>* node,
        const std::vector<std::string>& names)
    {
        int success = 1;
        for (size_t i = 0; i < names.size(); i++)
        {
            success = (success && (node->first_node(names[i].c_str()) != 0
                || node->first_attribute(names[i].c_str()) != 0));
        }
        return success;
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Check if values of a group of nodes or attributes are numbers.
    ////////////////////////////////////////////////////////////////////////////////
    int ReadIn::check_is_number(const rapidxml::xml_node<>* node,
        const std::vector<std::string>& names)
    {
        int success = 1;
        for (size_t i = 0; i < names.size(); i++)
        {
            if (node->first_node(names[i].c_str()))
                success = (success && (is_number(node->first_node(names[i].c_str())->value())));
            else if (node->first_attribute(names[i].c_str()))
                success = success && (is_number(node->first_attribute(names[i].c_str())->value()));
        }
        return success;
    }
    ////////////////////////////////////////////////////////////////////////////////
    // Check if inputs are numbers.
    ////////////////////////////////////////////////////////////////////////////////
    bool ReadIn::is_number(const std::string& token)
    {
        return std::regex_match(token, std::regex("[+\\-]?(?:0|[1-9]\\d*)(?:\\.\\d*)?(?:[eE][+\\-]?\\d+)?$"));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Check if inputs are integers.
    ////////////////////////////////////////////////////////////////////////////////
    bool ReadIn::is_integer(const std::string& token)
    {
        return std::regex_match(token, std::regex("(?:0|[1-9]\\d*)?$"));
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Check whether each of the string in a whitespace separated text is a number.
    ////////////////////////////////////////////////////////////////////////////////
    int ReadIn::extract_multi_number(const std::string token,
        std::vector<std::string>& number_string)
    {
        std::stringstream line_stream;
        line_stream << token;
        number_string.clear();
        string cell;
        int success = 1;
        while (line_stream >> cell)
        {
            number_string.push_back(cell);
            success = success && is_number(cell);
        }
        return success;
    }

} // End of namespace rope.
