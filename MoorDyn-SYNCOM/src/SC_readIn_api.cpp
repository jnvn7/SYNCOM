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


#include "SC_readIn_api.h"

using namespace std;
using namespace rapidxml;

namespace rope {

    ////////////////////////////////////////////////////////////////////////////////
    /// Initialization of the materDef.
    ////////////////////////////////////////////////////////////////////////////////
    ErrorCode ReadIn::readIn_data(string line_type, stressSolver& stressSolver)
    {
        // Get the absolute path to *input.xml file and set Output and Log files' location;
        std::string file_name = stressSolver.sc_folder + stressSolver.sc_inputfile;
        stressSolver.log_filename = stressSolver.sc_folder + "SynCOM_Log.txt";
        
        ////////////////////////////////////////////////////////////////////////////
        // MaterDef file.
        ////////////////////////////////////////////////////////////////////////////
        int flag = check_file_existence(file_name); 
        if (flag) 
            return ErrorCode::MATERDEF_FILE_NONEXISTENT;
        
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
            return ErrorCode::MATERDEF_FILE_ERROR_PARSE;
        }

        xml_document<> set_doc;
        std::string content(buffer.str());
        set_doc.parse<0>(&content[0]);

        xml_node<>* root_node;
        root_node = set_doc.first_node("MaterDef");
        if (root_node == 0) 
            return ErrorCode::MATERDEF_FILE_NO_MATERDEF_NODE; 
        
        ////////////////////////////////////////////////////////////////////////////
        // Check Material Properties Input Data.
        ////////////////////////////////////////////////////////////////////////////
        // Material properties;
        xml_node<>* child_node = root_node->first_node(line_type.c_str());
        if (child_node == 0)
            return ErrorCode::MATERDEF_FILE_MATERIAL_NAME_PROPERTIES;

        std::vector<std::string> names;
        names.push_back("sigma_yield0");
        names.push_back("MBL");
        names.push_back("Do");
        names.push_back("Dn");

        if (!check_availability(child_node, names))
            return ErrorCode::MATERDEF_FILE_INCOMPLETE_MATERIAL_PROPERTIES;
        else
        {
            names.erase(names.end() - 1);
            if (!check_is_number(child_node, names))
                return ErrorCode::MATERDEF_FILE_NAN_MATERIAL_PROPERTIES;
        }

        names.clear();

        stressSolver.material_props->sigma_yield0 =
            stod(child_node->first_node("sigma_yield0")->value());

        stressSolver.material_props->MBL =
            stod(child_node->first_node("MBL")->value());

        stressSolver.material_props->Do =
            stod(child_node->first_node("Do")->value());

        /// Read in Dn's values;
        std::vector<std::string> DnString;
        flag = extract_vector_element(child_node->first_node("Dn")->value(), DnString);

        if (DnString.size() <= 0)
            return ErrorCode::MATERDEF_FILE_BAD_DN_VALUES;
        else
        {
            stressSolver.material_props->Dn.resize(DnString.size());
            stressSolver.material_props->lamdaN.resize(DnString.size());
            for (int i = 0; i < DnString.size(); i++) {
                if (!is_number(DnString[i]))
                    return ErrorCode::NAN_INPUT_DATA;
                else {
                    stressSolver.material_props->Dn[i] = stod(DnString[i]);
                    if (stod(DnString[i]) < 0)
                        return ErrorCode::BAD_MATERIAL_PROPERTIES_INPUT;

                    stressSolver.material_props->lamdaN[i] = pow(10, -i);
                    stressSolver.material_props->sumDn += stressSolver.material_props->Dn[i];
                }
            }
        }
        
        // Clear parameters;
        DnString.clear();
        
        /// Read in a0, g0, g1, g2, Ep, np, H's values;
        std::vector<std::string> a0String;
        std::vector<std::string> g0String;
        std::vector<std::string> g1String;
        std::vector<std::string> g2String;
        std::vector<std::string> EpString;
        std::vector<std::string> npString;
        std::vector<std::string> HString;

        std::vector<std::string> a0stress_lim;
        std::vector<std::string> g0stress_lim;
        std::vector<std::string> g1stress_lim;
        std::vector<std::string> g2stress_lim;
        std::vector<std::string> Epstress_lim;
        std::vector<std::string> npstress_lim;
        std::vector<std::string> Hstress_lim;
        
        extract_vector_element(child_node->first_node("a0")->value(), a0String, a0stress_lim, stressSolver.material_props->step_num[0]);
        extract_vector_element(child_node->first_node("g0")->value(), g0String, g0stress_lim, stressSolver.material_props->step_num[1]);
        extract_vector_element(child_node->first_node("g1")->value(), g1String, g1stress_lim, stressSolver.material_props->step_num[2]);
        extract_vector_element(child_node->first_node("g2")->value(), g2String, g2stress_lim, stressSolver.material_props->step_num[3]);
        extract_vector_element(child_node->first_node("Ep")->value(), EpString, Epstress_lim, stressSolver.material_props->step_num[4]);
        extract_vector_element(child_node->first_node("np")->value(), npString, npstress_lim, stressSolver.material_props->step_num[5]);
        extract_vector_element(child_node->first_node("H")->value(), HString, Hstress_lim, stressSolver.material_props->step_num[6]);

        stressSolver.material_props->a0stress_lim.resize(2, std::vector<double>(a0stress_lim.size() / 2));
        stressSolver.material_props->g0stress_lim.resize(2, std::vector<double>(g0stress_lim.size() / 2));
        stressSolver.material_props->g1stress_lim.resize(2, std::vector<double>(g1stress_lim.size() / 2));
        stressSolver.material_props->g2stress_lim.resize(2, std::vector<double>(g2stress_lim.size() / 2));
        stressSolver.material_props->Epstress_lim.resize(2, std::vector<double>(Epstress_lim.size() / 2));
        stressSolver.material_props->npstress_lim.resize(2, std::vector<double>(npstress_lim.size() / 2));
        stressSolver.material_props->Hstress_lim.resize(2, std::vector<double>(Hstress_lim.size() / 2));
        
        if (a0String.size() <= 0 || g0String.size() <= 0 || g1String.size() <= 0
            || g2String.size() <= 0 || EpString.size() <= 0 || npString.size() <= 0
            || HString.size() <= 0)
            return ErrorCode::MATERDEF_FILE_NAN_MATERIAL_PROPERTIES;
        else
        {
            stressSolver.material_props->a0Coefs.resize(a0String.size());
            stressSolver.material_props->g0Coefs.resize(g0String.size());
            stressSolver.material_props->g1Coefs.resize(g1String.size());
            stressSolver.material_props->g2Coefs.resize(g2String.size());
            stressSolver.material_props->EpCoefs.resize(EpString.size());
            stressSolver.material_props->npCoefs.resize(npString.size());
            stressSolver.material_props->H_vpCoefs.resize(HString.size());
            
            // a0;
            for (size_t i = 0; i < a0String.size(); i++) {
                if (!is_number(a0String[i]))
                    return ErrorCode::NAN_INPUT_DATA;
                else
                    stressSolver.material_props->a0Coefs[i] = stod(a0String[i]);
            }
            for (size_t i = 0; i < a0stress_lim.size() / 2; i++) {
                if (!is_number(a0stress_lim[2 * i]) || !is_number(a0stress_lim[2 * i + 1]))
                    return ErrorCode::NAN_INPUT_DATA;
                else {
                    stressSolver.material_props->a0stress_lim[0][i] = stod(a0stress_lim[2 * i]);
                    stressSolver.material_props->a0stress_lim[1][i] = stod(a0stress_lim[2 * i + 1]);
                }
            }
            
            // g0;
            for (size_t i = 0; i < g0String.size(); i++) {
                if (!is_number(g0String[i]))
                    return ErrorCode::NAN_INPUT_DATA;
                else
                    stressSolver.material_props->g0Coefs[i] = stod(g0String[i]);
            }
            for (size_t i = 0; i < g0stress_lim.size() / 2; i++) {
                if (!is_number(g0stress_lim[2 * i]) || !is_number(g0stress_lim[2 * i + 1]))
                    return ErrorCode::NAN_INPUT_DATA;
                else {
                    stressSolver.material_props->g0stress_lim[0][i] = stod(g0stress_lim[2 * i]);
                    stressSolver.material_props->g0stress_lim[1][i] = stod(g0stress_lim[2 * i + 1]);
                }
            }

            // g1;
            for (size_t i = 0; i < g1String.size(); i++) {
                if (!is_number(g1String[i]))
                    return ErrorCode::NAN_INPUT_DATA;
                else
                    stressSolver.material_props->g1Coefs[i] = stod(g1String[i]);
            }
            for (size_t i = 0; i < g1stress_lim.size() / 2; i++) {
                if (!is_number(g1stress_lim[2 * i]) || !is_number(g1stress_lim[2 * i + 1]))
                    return ErrorCode::NAN_INPUT_DATA;
                else {
                    stressSolver.material_props->g1stress_lim[0][i] = stod(g1stress_lim[2 * i]);
                    stressSolver.material_props->g1stress_lim[1][i] = stod(g1stress_lim[2 * i + 1]);
                }
            }

            // g2;
            for (size_t i = 0; i < g2String.size(); i++) {
                if (!is_number(g2String[i]))
                    return ErrorCode::NAN_INPUT_DATA;
                else
                    stressSolver.material_props->g2Coefs[i] = stod(g2String[i]);
            }
            for (size_t i = 0; i < g2stress_lim.size() / 2; i++) {
                if (!is_number(g2stress_lim[2 * i]) || !is_number(g2stress_lim[2 * i + 1]))
                    return ErrorCode::NAN_INPUT_DATA;
                else {
                    stressSolver.material_props->g2stress_lim[0][i] = stod(g2stress_lim[2 * i]);
                    stressSolver.material_props->g2stress_lim[1][i] = stod(g2stress_lim[2 * i + 1]);
                }
            }

            // Ep;
            for (size_t i = 0; i < EpString.size(); i++) {
                if (!is_number(EpString[i]))
                    return ErrorCode::NAN_INPUT_DATA;
                else
                    stressSolver.material_props->EpCoefs[i] = stod(EpString[i]);
            }
            for (size_t i = 0; i < Epstress_lim.size() / 2; i++) {
                if (!is_number(Epstress_lim[2 * i]) || !is_number(Epstress_lim[2 * i + 1]))
                    return ErrorCode::NAN_INPUT_DATA;
                else {
                    stressSolver.material_props->Epstress_lim[0][i] = stod(Epstress_lim[2 * i]);
                    stressSolver.material_props->Epstress_lim[1][i] = stod(Epstress_lim[2 * i + 1]);
                }
            }

            // np;
            for (size_t i = 0; i < npString.size(); i++) {
                if (!is_number(npString[i]))
                    return ErrorCode::NAN_INPUT_DATA;
                else
                    stressSolver.material_props->npCoefs[i] = stod(npString[i]);
            }
            for (size_t i = 0; i < npstress_lim.size() / 2; i++) {
                if (!is_number(npstress_lim[2 * i]) || !is_number(npstress_lim[2 * i + 1]))
                    return ErrorCode::NAN_INPUT_DATA;
                else {
                    stressSolver.material_props->npstress_lim[0][i] = stod(npstress_lim[2 * i]);
                    stressSolver.material_props->npstress_lim[1][i] = stod(npstress_lim[2 * i + 1]);
                }
            }

            // H_vp;
            for (size_t i = 0; i < HString.size(); i++) {
                if (!is_number(HString[i]))
                    return ErrorCode::NAN_INPUT_DATA;
                else
                    stressSolver.material_props->H_vpCoefs[i] = stod(HString[i]);
            }
            for (size_t i = 0; i < Hstress_lim.size() / 2; i++) {
                if (!is_number(Hstress_lim[2 * i]) || !is_number(Hstress_lim[2 * i + 1]))
                    return ErrorCode::NAN_INPUT_DATA;
                else {
                    stressSolver.material_props->Hstress_lim[0][i] = stod(Hstress_lim[2 * i]);
                    stressSolver.material_props->Hstress_lim[1][i] = stod(Hstress_lim[2 * i + 1]);
                }
            }
           
            // Clear parameters;
            a0String.clear(); a0stress_lim.clear(); g0String.clear(); g0stress_lim.clear();
            g1String.clear(); g1stress_lim.clear(); g2String.clear(); g2stress_lim.clear();
            EpString.clear(); Epstress_lim.clear(); npString.clear(); npstress_lim.clear();
            HString.clear(); Hstress_lim.clear();
        }
    
        ////////////////////////////////////////////////////////////////////////////
        // Check Material Properties Input Data.
        ////////////////////////////////////////////////////////////////////////////
        child_node = root_node->first_node("numerical_setting");
        if (child_node == 0)
            return ErrorCode::MATERDEF_FILE_MATERIAL_NAME_PROPERTIES;

        names.push_back("limit");
        names.push_back("tol");

        if (!check_availability(child_node, names))
            return ErrorCode::MATERDEF_FILE_INCOMPLETE_MATERIAL_PROPERTIES;
        else
        {
            if (!check_is_number(child_node, names))
                return ErrorCode::MATERDEF_FILE_NAN_MATERIAL_PROPERTIES;
        }

        stressSolver.material_props->limit = stoi(child_node->first_node("limit")->value());

        stressSolver.material_props->tol = stod(child_node->first_node("tol")->value());

        // Clear parameters;
        names.clear();

        return ErrorCode::SUCCESS;
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
    /// Extract elements in input vector(string).
    ////////////////////////////////////////////////////////////////////////////////
    int ReadIn::extract_vector_element(const std::string token,
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

    ////////////////////////////////////////////////////////////////////////////////
    /// Extract polynomial coefficients (with step limits) from the input strings.
    ////////////////////////////////////////////////////////////////////////////////
    int ReadIn::extract_vector_element(const std::string token,
        std::vector<std::string>& number_string, std::vector<std::string>& step_lim,
        int& step_num)
    {
        std::stringstream line_stream;
        line_stream << token;
        number_string.clear();
        step_lim.clear();
        string cell;
        char limit = 'L';
        int success = 1, count = 0;
        while (line_stream >> cell)
        {
            if (limit != cell.at(0)) {
                number_string.push_back(cell);
                success = success && is_number(cell);
                count += 1;
            }
            else {
                cell.erase(0, 2); cell.erase(cell.end() - 1);
                step_lim.push_back(cell);
                step_lim.push_back(to_string(count));
                step_num += 1;
            }
        }
        return success;
    }

} // End of namespace rope.
