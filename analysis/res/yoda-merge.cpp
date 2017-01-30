#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
namespace po = boost::program_options;
using namespace boost;
using namespace std;

const char_separator<char> separator{"\t"};

enum class HistType { Histo1D, Histo2D };

struct histRecord {
  HistType hist_type;
  double xlow = 0.;
  double xhigh = 0.;
  double ylow = 0.;
  double yhigh = 0.;
  double sumw = 0.;
  double sumw2 = 0.;
  double sumwx = 0.;
  double sumwx2 = 0.;
  double sumwy = 0.;
  double sumwy2 = 0.;
  double sumwxy = 0.;
  unsigned long numEntries = 0;
};

void readRecord(histRecord &rec, const string &line) {
  tokenizer<char_separator<char>> tok(line, separator);
  auto iter = tok.begin();
  if (rec.hist_type == HistType::Histo1D) {
    rec.xlow = stod(*iter++);
    rec.xhigh = stod(*iter++);
    rec.sumw += stod(*iter++);
    rec.sumw2 += stod(*iter++);
    rec.sumwx += stod(*iter++);
    rec.sumwx2 += stod(*iter++);
    rec.numEntries += stoul(*iter++);
  } else if (rec.hist_type == HistType::Histo2D) {
    rec.xlow = stod(*iter++);
    rec.xhigh = stod(*iter++);
    rec.ylow = stod(*iter++);
    rec.yhigh = stod(*iter++);
    rec.sumw += stod(*iter++);
    rec.sumw2 += stod(*iter++);
    rec.sumwx += stod(*iter++);
    rec.sumwx2 += stod(*iter++);
    rec.sumwy += stod(*iter++);
    rec.sumwy2 += stod(*iter++);
    rec.sumwxy += stod(*iter++);
    rec.numEntries += stoul(*iter++);
  } else {
    cerr << "COULD NOT PARSE LINE: " << line << "\n";
    throw std::runtime_error("COULD NOT PARSE LINE");
  }
}

struct statsRecord {
  double sumw = 0.;
  double sumw2 = 0.;
  double sumwx = 0.;
  double sumwx2 = 0.;
  double sumwy = 0.;
  double sumwy2 = 0.;
  double sumwxy = 0.;
  unsigned long numEntries = 0;

  statsRecord &operator+=(const statsRecord &other);
};

statsRecord &statsRecord::operator+=(const statsRecord &other) {
  sumw += other.sumw;
  sumw2 += other.sumw2;
  sumwx += other.sumwx;
  sumwx2 += other.sumwx2;
  sumwy += other.sumwy;
  sumwy2 += other.sumwy2;
  sumwxy += other.sumwxy;
  numEntries += other.numEntries;

  return *this;
}

struct header {
  HistType hist_type;
  string path;
  string title;
  statsRecord total;
  statsRecord underflow;
  statsRecord overflow;

  header(ifstream &input, int fileNum);
  header() = default;
  header(const header &) = default;
  header &operator=(const header &) = default;
  header &operator+=(const header &other);
  bool operator!=(const header &other);
};

header &header::operator+=(const header &other) {
  total += other.total;
  underflow += other.underflow;
  overflow += other.overflow;
  return *this;
}

bool header::operator!=(const header &other) {
  return !((path == other.path) && (title == other.title));
}

header::header(ifstream &input, int fileNum) {
  if (!input.good())
    throw std::runtime_error("Error reading input file");
  if (input.tellg() != 0) {
    cerr << "Warning: reading header but not at start of file.\n"
         << "Returning to start of file.\n";
    input.seekg(0);
  }

  string buffer;
  if (!getline(input, buffer))
    throw std::runtime_error("File Empty");
  if (buffer.substr(0, 18) == "BEGIN YODA_HISTO1D")
    hist_type = HistType::Histo1D;
  else if (buffer.substr(0, 18) == "BEGIN YODA_HISTO2D")
    hist_type = HistType::Histo2D;
  else
    throw std::runtime_error("Not a YODA histogram");
  // A lot less resilient from here
  getline(input, buffer);
  path = buffer.substr(5, 1005); // path is presumably < 1000 chars
  getline(input, buffer);
  title = buffer.substr(6, 1006); // title is presumably < 1000 chars
  getline(input, buffer);         // skip type line
  // Total (skip comments)
  getline(input, buffer);
  while (buffer.substr(0, 1) == "#")
    getline(input, buffer);
  try {
    tokenizer<char_separator<char>> total_tok(buffer, separator);
    int total_idx = 0;
    for (const auto &token : total_tok) {
      try {
        switch (total_idx) {
        case 2:
          total.sumw = stod(token);
          break;
        case 3:
          total.sumw2 = stod(token);
          break;
        case 4:
          total.sumwx = stod(token);
          break;
        case 5:
          total.sumwx2 = stod(token);
          break;
        case 6:
          if (hist_type == HistType::Histo1D)
            total.numEntries = stoul(token);
          else if (hist_type == HistType::Histo2D)
            total.sumwy = stod(token);
          break;
        // below cases only appear for Histo2D
        case 7:
          total.sumwy2 = stod(token);
          break;
        case 8:
          total.sumwxy = stod(token);
          break;
        case 9:
          total.numEntries = stoul(token);
          break;
        }
      } catch (...) {
        cerr << "TOKEN WAS " << token << "\n";
        throw;
      }
      ++total_idx;
    }
  } catch (...) {
    cerr << "FILE " << fileNum;
    throw std::runtime_error("EXCEPTION PARSING HEADER TOTAL\n");
  }
  // Underflow
  if (hist_type == HistType::Histo1D) {
    getline(input, buffer);
    try {
      tokenizer<char_separator<char>> underflow_tok(buffer, separator);
      int underflow_idx = 0;
      for (const auto &token : underflow_tok) {
        switch (underflow_idx) {
        case 2:
          underflow.sumw = stod(token);
          break;
        case 3:
          underflow.sumw2 = stod(token);
          break;
        case 4:
          underflow.sumwx = stod(token);
          break;
        case 5:
          underflow.sumwx2 = stod(token);
          break;
        case 6:
          underflow.numEntries = stoul(token);
          break;
        }
        ++underflow_idx;
      }
    } catch (...) {
      cerr << "FILE " << fileNum;
      throw std::runtime_error("EXCEPTION PARSING HEADER UNDERFLOW\n");
    }
    // Overflow
    getline(input, buffer);
    try {
      tokenizer<char_separator<char>> overflow_tok(buffer, separator);
      int overflow_idx = 0;
      for (const auto &token : overflow_tok) {
        switch (overflow_idx) {
        case 2:
          overflow.sumw = stod(token);
          break;
        case 3:
          overflow.sumw2 = stod(token);
          break;
        case 4:
          overflow.sumwx = stod(token);
          break;
        case 5:
          overflow.sumwx2 = stod(token);
          break;
        case 6:
          overflow.numEntries = stoul(token);
          break;
        }
        ++overflow_idx;
      }
    } catch (...) {
      cerr << "FILE " << fileNum;
      throw std::runtime_error("EXCEPTION PARSING HEADER OVERFLOW\n");
    }
  }
  // skip line
  getline(input, buffer);
  // and another for Histo2D
  if (hist_type == HistType::Histo2D)
    getline(input, buffer);
}

ostream &operator<<(ostream &os, const header &head) {
  auto &&total = head.total;
  auto &&underflow = head.underflow;
  auto &&overflow = head.overflow;
  os << std::scientific << std::left;
  if (head.hist_type == HistType::Histo1D)
    os << "BEGIN YODA_HISTO1D ";
  else if (head.hist_type == HistType::Histo2D)
    os << "BEGIN YODA_HISTO2D ";
  os << head.path << "\n";
  os << "Path=" << head.path << "\n";
  os << "Title=" << head.title << "\n";
  os << "Type=" << (head.hist_type == HistType::Histo1D ? "Histo1D" : "Histo2D")
     << "\n"
     << "# elided\n"
     << (head.hist_type == HistType::Histo1D ? "# elided\n" : "");
  // Histo1D has 2 comments here ^^^
  os << "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\t";
  if (head.hist_type == HistType::Histo2D)
    os << "sumwy\tsumwy2\tsumwxy\t";
  os << "numEntries\n";
  os << setw(8) << "Total"
     << "\t" << setw(8) << "Total"
     << "\t" << total.sumw << "\t" << total.sumw2 << "\t" << total.sumwx << "\t"
     << total.sumwx2 << "\t";
  if (head.hist_type == HistType::Histo2D)
    os << total.sumwy << "\t" << total.sumwy2 << "\t" << total.sumwxy << "\t";
  os << total.numEntries << "\n";
  if (head.hist_type == HistType::Histo1D) {
    os << "Underflow\tUnderflow\t" << underflow.sumw << "\t" << underflow.sumw2
       << "\t" << underflow.sumwx << "\t" << underflow.sumwx2 << "\t"
       << underflow.numEntries << "\n";
    os << "Overflow\tOverflow\t" << overflow.sumw << "\t" << overflow.sumw2
       << "\t" << overflow.sumwx << "\t" << overflow.sumwx2 << "\t"
       << overflow.numEntries << "\n";
    os << "# xlow\txhigh\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries"
       << "\n";
  } else if (head.hist_type == HistType::Histo2D) {
    os << "# 2D outflow persistency not currently supported until API is "
          "stable\n";
    os << "# "
          "xlow\txhigh\tylow\tyhigh\tsumw\tsumw2\tsumwx\tsumwx2\tsumwy\tsumwy2"
          "\tsumwxy\tnumEntries\n";
  }
  return os;
}

ostream &operator<<(ostream &os, const histRecord &rec) {
  os << std::scientific;
  os << rec.xlow << "\t" << rec.xhigh << "\t";
  if (rec.hist_type == HistType::Histo2D)
    os << rec.ylow << "\t" << rec.yhigh << "\t";
  os << rec.sumw << "\t" << rec.sumw2 << "\t";
  os << rec.sumwx << "\t" << rec.sumwx2 << "\t";
  if (rec.hist_type == HistType::Histo2D)
    os << rec.sumwy << "\t" << rec.sumwy2 << "\t" << rec.sumwxy << "\t";
  os << rec.numEntries;

  return os;
}

int main(int argc, char **argv) {
  string outputFilename{};
  string flatOutputFilename{};
  vector<string> inputFilenames{};

  po::options_description opts("Merge YODA histograms");
  opts.add_options()("help,h", "print help")(
      "output,o", po::value<string>(&outputFilename), "YODA output file")(
      "flat-output,f", po::value<string>(&flatOutputFilename),
      "FLAT output file")("input,i", po::value<vector<string>>(&inputFilenames),
                          "input files");
  po::positional_options_description posOpts;
  posOpts.add("input", -1);
  po::variables_map varMap;
  po::store(po::command_line_parser(argc, argv)
                .options(opts)
                .positional(posOpts)
                .run(),
            varMap);
  po::notify(varMap);
  if (varMap.count("help")) {
    // print help
    cout << opts << endl;
    return 0;
  } else if (!varMap.count("input")) {
    cerr << opts << "\n";
    cerr << "Error: no input files listed" << endl;
    return 1;
  } else if (!varMap.count("output")) {
    cerr << opts << "\n";
    cerr << "Error: no output files listed" << endl;
    return 1;
  }

  // Start combining files
  vector<ifstream> inputFiles;
  for (const auto &file : inputFilenames) {
    inputFiles.emplace_back(file);
    if (!inputFiles.back().good()) {
      cerr << "Failed to open " << file << endl;
      return 2;
    }
  }

  ofstream outputFile{outputFilename};
  if (!outputFile.good()) {
    cerr << "Failed to open output file" << endl;
    return 2;
  }

  header outputHeader;
  vector<histRecord> histogram{};
  bool firstFile = true; // false after processing first file
  int fileNum = -1;
  for (auto &&file : inputFiles) {
    ++fileNum;
    if (firstFile) {
      outputHeader = header(file, fileNum);
      int recIdx = -1;
      for (string buffer;
           getline(file, buffer), buffer.substr(0, 3) != "END";) {
        ++recIdx;
        histogram.emplace_back();
        histRecord &rec = histogram.back();
        rec.hist_type = outputHeader.hist_type;
        try {
          readRecord(rec, buffer);
        } catch (const std::exception &except) {
          cerr << "EXCEPTION: " << except.what() << "\n";
          cerr << "WHILE PROCESSING FILE " << fileNum << " RECORD " << recIdx
               << endl;
          return 5;
        }
      }
      firstFile = false;
    } else {
      // ALL BUT FIRST FILE
      header head(file, fileNum);
      if (head != outputHeader) {
        cerr << "Error: names / titles do not all match (file " << fileNum
             << ")" << endl;
        return 3;
      }
      outputHeader += head;
      auto iter = histogram.begin();
      int recIdx = -1;
      for (string buffer;
           getline(file, buffer), buffer.substr(0, 3) != "END";) {
        ++recIdx;
        if (iter == histogram.end()) {
          cerr << "Error: histogram lengths do not all match (file " << fileNum
               << ")" << endl;
          return 3;
        }
        histRecord &rec = *iter;
        rec.hist_type = head.hist_type;
        try {
          readRecord(rec, buffer);
        } catch (const std::exception except) {
          cerr << "EXCEPTION: " << except.what() << "\n";
          cerr << "WHILE PROCESSING FILE " << fileNum << "RECORD " << recIdx
               << endl;
          return 5;
        }
        ++iter;
      }
    }
  }
  // Write out merged histogram

  outputFile << outputHeader;
  for (auto &&rec : histogram) {
    outputFile << rec << "\n";
  }
  if (outputHeader.hist_type == HistType::Histo1D)
    outputFile << "END YODA_HISTO1D\n\n" << flush;
  else if (outputHeader.hist_type == HistType::Histo2D)
    outputFile << "END YODA_HISTO2D\n\n" << flush;

  // write out flat file
  if (varMap.count("flat-output")) {
    ofstream os{flatOutputFilename};
    if (!os.good()) {
      cerr << "Failed to open flat output file";
      return 10;
    }
    if (outputHeader.hist_type == HistType::Histo1D) {
      os << "# BEGIN HISTO1D " << outputHeader.path << "\n";
      os << "Path=" << outputHeader.path << "\n";
      os << "Title=" << outputHeader.title << "\n";
      os << "# xlow\txhigh\tval\terrminus\terrplus\n";

      for (auto &&rec : histogram) {
        os << std::scientific;
        os << rec.xlow << "\t" << rec.xhigh << "\t";
        os << (rec.sumw / (rec.xhigh - rec.xlow)) << "\t";
        double error = sqrt(rec.sumw2) / (rec.xhigh - rec.xlow);
        os << error << "\t" << error << "\n";
      }
      os << "# END HISTO1D\n\n";
    } else if (outputHeader.hist_type == HistType::Histo2D) {
      os << "# BEGIN HISTO2D " << outputHeader.path << "\n";
      os << "Path=" << outputHeader.path << "\n";
      os << "Title=" << outputHeader.title << "\n";
      os << "# xlow\txhigh\tylow\tyhigh\tval\terrminus\terrplus\n";

      for (auto &&rec : histogram) {
        os << std::scientific;
        os << rec.xlow << "\t" << rec.xhigh << "\t" << rec.ylow << "\t"
           << rec.yhigh << "\t";
        os << (rec.sumw / ((rec.xhigh - rec.xlow) * (rec.yhigh - rec.ylow)))
           << "\t";
        double error =
            sqrt(rec.sumw2) / ((rec.xhigh - rec.xlow) * (rec.yhigh - rec.ylow));
        os << error << "\t" << error << "\n";
      }
      os << "# END HISTO2D\n\n";
    }
  }
  return 0;
}
