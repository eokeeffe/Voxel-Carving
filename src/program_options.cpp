/*
   This example of boost program options parser only analyzes command line options

   g++ program_options.cpp -l boost_program_options

 */
#include <iostream>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
    // Declare the supported options.
    int var2;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("var1", po::value<int>(), "set var1")
        ("var2", po::value<int>(&var2)->default_value(10), "set var2")
        ("string1", po::value<string>(), "set string1")
        ("list1", po::value< vector<string> >(), "set list1")
        ("input-file", po::value< vector<string> >(), "set input files")
        ;

    // Declare which options are positional
    po::positional_options_description p;
    p.add("input-file", -1);

    po::variables_map vm;
    // parse regular options
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    // parse positional options
    po::store(po::command_line_parser(argc, argv). options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << endl;
        return 1;
    }

    if (vm.count("var1")) {
        cout << "var1 was set to " 
            << vm["var1"].as<int>() << ".\n";
    } else {
        cout << "var1 was not set.\n";
    }

    if (vm.count("var2")) {
        cout << "var2 was set to " 
            << vm["var2"].as<int>() << ".\n";
    }

    if (vm.count("string1")) {
        cout << "string1 was set to " << vm["string1"].as<string>() << ".\n";
    }

    if (vm.count("list1")) {
        cout << "list1 is: " << endl;
        for(int i = 0; i < vm["list1"].as< vector<string> >().size(); i++)
            cout << vm["list1"].as< vector<string> >()[i] << endl;
    }

    if (vm.count("input-file")) {
        cout << "input files are: " << endl;
        for(int i = 0; i < vm["input-file"].as< vector<string> >().size(); i++)
            cout << vm["input-file"].as< vector<string> >()[i] << endl;
    }


    return EXIT_SUCCESS;
}