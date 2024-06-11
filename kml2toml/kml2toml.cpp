#include <tinyxml2.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <string>

void process_whitnesses(tinyxml2::XMLElement*, std::fstream&);
void process_placemark(tinyxml2::XMLElement*, int, std::fstream&);
const char* last_word(const std::string&);
const char* last_word_dur(const std::string&);

int main(int argc, char** argv) {
    // Check args
    if (argc < 2) {
        std::cerr << "Usage: kml file1 [file2 ...]" << std::endl;
        return 1;
    }

    // Open output file
    std::fstream outfile;
    outfile.open("kml.toml", std::ios::out);
    if (!outfile) {
        std::cerr << "Failed to open 'kml.toml'";
        return 1;
    }

    // Run for each file
    for (int i = 1; i < argc; i++) {
        // Open document
        tinyxml2::XMLDocument doc(true, tinyxml2::COLLAPSE_WHITESPACE);
        doc.LoadFile(argv[i]);

        // Let's make sure the file loaded fine...
        if (doc.ErrorID()) {
            std::cerr << "Error: Failed to open file\n";
            break;
        }

        // For each "Folder"
        tinyxml2::XMLElement* pRoot =
            doc.RootElement()->FirstChildElement()->FirstChildElement("Folder");
        while (pRoot) {
            const char* folder = pRoot->FirstChildElement("name")->GetText();
            if (!strcmp(folder, "Witnesses"))
                process_whitnesses(pRoot, outfile);
            pRoot = pRoot->NextSiblingElement("Folder");
        }
    }

    // Answer template (filled manually)
    outfile << "[answer]\n";
    outfile << "lat = \n";
    outfile << "lon = \n";
    outfile << "h = \n";
    outfile << "bearing = \n";
    outfile << "incidence = \n";
    outfile << "speed = \n";
    outfile << "\n";
    outfile << "[meta]\n";
    outfile << "url_imo = \"\"\n";
    outfile << "url_strewnify = \"\"\n";

    // Close file and exit
    outfile.close();
    return 0;
}

void process_whitnesses(tinyxml2::XMLElement* pRoot, std::fstream& outfile) {
    pRoot = pRoot->FirstChildElement("Folder");
    while (pRoot) {
        // Folder name is 'Experience Level x'
        const char* folder = pRoot->FirstChildElement("name")->GetText();
        int experience = atoi(folder + 17);

        // Go through 'Placemark'
        tinyxml2::XMLElement* pPlacemark = pRoot->FirstChildElement("Placemark");
        while (pPlacemark) {
            process_placemark(pPlacemark, experience, outfile);
            pPlacemark = pPlacemark->NextSiblingElement("Placemark");
        }

        // Next folder
        pRoot = pRoot->NextSiblingElement("Folder");
    }
}

void process_placemark(tinyxml2::XMLElement* pPlacemark, int experience, std::fstream& outfile) {
    // Get basic info
    const char* name = pPlacemark->FirstChildElement("name")->GetText();
    std::string desc(pPlacemark->FirstChildElement("description")->GetText());
    double lon, lat;
    sscanf(
        pPlacemark->FirstChildElement("Point")->FirstChildElement("coordinates")->GetText(),
        "%lf,%lf", &lon, &lat
    );

    // Height
    std::regex height_regex("([0-9]+\\.)?[0-9]+m");
    std::smatch height_match;
    double h = 0;
    if (!std::regex_search(desc, height_match, height_regex)) {
        std::cout << "Observer '" << name << "' does not have 'height', using '0'\n";
        std::cout << "\tLocation: " << lat << "," << lon << "\n";
    } else
        h = std::stof(height_match[0]);

    // Desent Angle
    std::regex da_regex("Descent Angle </th> <td[^<>]*>\\W*[0-9\\.]+");
    std::smatch da_match;
    double a = -1;
    if (std::regex_search(desc, da_match, da_regex))
        sscanf(last_word(da_match[0].str()), "%lf", &a);

    // Azimuth Begin
    std::regex zb_regex("First azimuth </th> <td[^<>]*>\\W*[0-9\\.]+");
    std::smatch zb_match;
    double zb = -1;
    if (std::regex_search(desc, zb_match, zb_regex))
        sscanf(last_word(zb_match[0].str()), "%lf", &zb);

    // Altitude Begin
    std::regex hb_regex("First elevation </th> <td[^<>]*>\\W*[0-9\\.]+");
    std::smatch hb_match;
    double hb = -1;
    if (std::regex_search(desc, hb_match, hb_regex))
        sscanf(last_word(hb_match[0].str()), "%lf", &hb);

    // Azimuth End
    std::regex z0_regex("Last azimuth </th> <td[^<>]*>\\W*[0-9\\.]+");
    std::smatch z0_match;
    double z0 = -1;
    if (std::regex_search(desc, z0_match, z0_regex))
        sscanf(last_word(z0_match[0].str()), "%lf", &z0);

    // Altitude End
    std::regex h0_regex("Last elevation </th> <td[^<>]*>\\W*[0-9\\.]+");
    std::smatch h0_match;
    double h0 = -1;
    if (std::regex_search(desc, h0_match, h0_regex))
        sscanf(last_word(h0_match[0].str()), "%lf", &h0);

    // Duration
    std::regex t_regex("Duration</th> </th> <td[^<>?]*>(&ap;|<)?[0-9\\.]+");
    std::smatch t_match;
    double t = -1;
    if (std::regex_search(desc, t_match, t_regex))
        sscanf(last_word_dur(t_match[0].str()), "%lf", &t);

    auto push_toml = [&](const char* label, double value) {
        if (value != -1)
            outfile << label << " = " << value << std::endl;
    };

    // Write to file
    outfile << "[[sample]]" << std::endl;
    outfile << "lat = " << lat << std::endl;
    outfile << "lon = " << lon << std::endl;
    outfile << "h = " << h << std::endl;
    push_toml("a", a);
    push_toml("zb", zb);
    push_toml("hb", hb);
    push_toml("z0", z0);
    push_toml("h0", h0);
    push_toml("t", t);
    outfile << "name = \"" << name << '"' << std::endl;
    outfile << "exp = " << experience << std::endl << std::endl;
}

const char* last_word(const std::string& str) {
    int i = str.size() - 1;
    while (i > 0 && str[i] != ' ' && str[i] != '\t')
        i--;
    return str.c_str() + i + 1;
}

const char* last_word_dur(const std::string& str) {
    int i = str.size() - 1;
    while (i > 0 && str[i] != ';' && str[i] != '<')
        i--;
    return str.c_str() + i + 1;
}
