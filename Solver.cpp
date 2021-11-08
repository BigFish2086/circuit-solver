#include "eigenLibrary/Eigen/Eigen"
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace Eigen;
using namespace std;

struct Component {
  string name;
  complex<double> voltage = 0;
  complex<double> current = 0;
  complex<double> power = 0;
  complex<double> Z = 0;
  complex<double> Zconj = 0;
  complex<double> Y = 0;
  // nodes entered as (-ve) for (node1) and (+ve) for (node2) ,espesially  for
  // VS and CS.
  int node1 = 0, node2 = 0, VS_num = 0;
  // for dependent sources;
  int posNode = 0, negNode = 0;
  complex<double> value = 0;
};

const double PI = 3.14159265359;
int w = 0; // frequency -> entered at the start of the text/input file.
// counters for circuit elements.
int numOfVS = 0, numOfCS = 0, numOfR = 0, numOfC = 0, numOfL = 0,
    numOfNodes = 0;
int numOfH = 0;      // number of C.C.V.S.
int numOfF = 0;      // number of C.C.C.S.
int numOfE = 0;      // number of V.C.V.S.
int numOfG = 0;      // number of V.C.C.S.
int numOfZeroVs = 0; // number of Zero Vs.
ifstream inputFile;
ofstream outputFile;
vector<Component> components;
vector<vector<Component>> nodesVector;
vector<vector<complex<double>>> nodePolarity; // positive or negative
vector<string> names;
Eigen::MatrixXcd A, G, B, C, D; // matrix A
Eigen::MatrixXcd Z, I, E;       // matrix Z
Eigen::MatrixXcd X;             // matrix X

void SplitBySpaces(const string &line, vector<string> &words) {
  string word = "";
  stringstream ss(line);
  while (getline(ss, word, ' ')) {
    words.push_back(word);
  }
}

bool IsOkParamteres(const vector<string> &words, const int &size) {
  if (words.size() == size)
    return true;
  else
    return false;
}

void FillSourcesData(const vector<string> &words, const int &type) {
  // type_0 is for VS and type_1 is for CS
  Component source;
  source.name = words[0];
  source.node1 = stoi(words[1]);
  source.node2 = stoi(words[2]);
  double magntiude = stod(words[3]);
  double phase = stod(words[4]);
  if (type == 0) {
    source.voltage = complex<double>(magntiude * cos(phase * PI / 180),
                                     magntiude * sin(phase * PI / 180));
    names.push_back(source.name);
  } else if (type == 1) {
    source.current = complex<double>(magntiude * cos(phase * PI / 180),
                                     magntiude * sin(phase * PI / 180));
  }
  components.push_back(source);
}

bool FillElementsData(const vector<string> &words, const int &type) {
  // type_0 is for R , type_1 is for C ,type_2 is for L.
  Component element;
  element.name = words[0];
  element.node1 = stoi(words[1]);
  element.node2 = stoi(words[2]);
  stringstream ss(words[3]);
  double magntiude = 0;
  string j = "";
  ss >> magntiude >> j;
  // cout << "mag = " << magntiude << " j = " << j << endl;
  magntiude = stod(words[3]);
  double infinity = 1000000;
  if (w > 0) {
    if (j == "") {
      if (type == 0) {
        element.Z = complex<double>(magntiude, 0);
      } else if (type == 1) {
        element.Z = complex<double>(0, -1 / (w * magntiude));
      } else if (type == 2) {
        element.Z = complex<double>(0, w * magntiude);
      }
    } else if (j == "j") {
      if (type == 0) {
        element.Z = complex<double>(magntiude, 0);
      } else if (type == 1) {
        element.Z = complex<double>(0, magntiude);
      } else if (type == 2) {
        element.Z = complex<double>(0, magntiude);
      }
    } else {
      cerr << "Element " << words[0] << " has wrong data." << endl;
      return false;
    }
    element.Y = pow(element.Z, -1);
    element.Zconj = conj(element.Z);
    components.push_back(element);
    return true;
  } else if (w == 0) {
    if (type == 0) {
      element.Z = complex<double>(magntiude, 0);
      element.Y = pow(element.Z, -1);
    } else if (type == 1) {
      element.Z = complex<double>(infinity, 0);
      element.Y = complex<double>(0, 0);
    } else if (type == 2) {
      element.Z = complex<double>(0, 0);
      element.Y = complex<double>(infinity, 0);
    }
    element.Zconj = conj(element.Z);
    components.push_back(element);
    return true;
  }
}

void FillDependantData(const vector<string> &words, const int &type) {
  // type_0 is for CCVS , type_1 is for CCCS.
  // type_2 is for VCVS , type_3 is for VCCS.
  Component depSource;
  depSource.name = words[0];
  if (type == 0 || type == 2) {
    names.push_back(depSource.name); //[TODO....]
  }
  depSource.node1 = stoi(words[1]);
  depSource.node2 = stoi(words[2]);
  depSource.negNode = stoi(words[3]);
  depSource.posNode = stoi(words[4]);
  depSource.value = stod(words[5]);
  components.push_back(depSource);
}

bool ReadCircuitData() {
  string fileName;
  cout << "Enter the file Name ";
  cin >> fileName;
  regex reg("([\w.-]*)(.txt)");
  if (!(std::regex_match(fileName, reg)))
    fileName += ".txt";
  inputFile.open(fileName);
  if (!inputFile.is_open()) {
    cerr << "Wrong file name." << endl;
    inputFile.close();
    return false;
  } else {
    string line;
    vector<string> words;
    getline(inputFile, line);
    SplitBySpaces(line, words);
    w = stod(words[0]);
    if (w < 0) {
      cerr << "omega w can't be less than zero.";
      return false;
    }
    while (!inputFile.eof()) {
      bool valid;
      words.clear();
      getline(inputFile, line);
      SplitBySpaces(line, words);
      if (words[0].find("Vs") != string::npos) { // normal voltage source
        valid = IsOkParamteres(words, 5);
        if (!valid) {
          cerr << "Voltage source has invalid parameters " << endl;
          inputFile.close();
          return false;
        } else {
          FillSourcesData(words, 0);
          numOfVS++;
        }
      } else if (words[0].find("Cs") != string::npos) { // normal current source
        valid = IsOkParamteres(words, 5);
        if (!valid) {
          cerr << "Current source has invalid parameters " << endl;
          inputFile.close();
          return false;
        } else {
          FillSourcesData(words, 1);
          numOfCS++;
        }
      } else if (words[0].find("R") != string::npos) { // Resistatnce
        valid = IsOkParamteres(words, 4);
        if (!valid) {
          cerr << "Resistor has invalid parameters " << endl;
          inputFile.close();
          return false;
        } else {
          if (!FillElementsData(words, 0)) {
            return false;
          }
          numOfR++;
        }
      } else if (words[0].find("C") != string::npos) { // Capacitor
        valid = IsOkParamteres(words, 4);
        if (!valid) {
          cerr << "Capacitor has invalid parameters " << endl;
          inputFile.close();
          return false;
        } else {
          if (!FillElementsData(words, 1)) {
            return false;
          }
          numOfC++;
        }
      } else if (words[0].find("L") != string::npos) { // Inductor
        valid = IsOkParamteres(words, 4);
        if (!valid) {
          cerr << "Inductor has invalid parameters " << endl;
          inputFile.close();
          return false;
        } else {
          if (!FillElementsData(words, 2)) {
            return false;
          }
          numOfL++;
        }
      } else if (words[0].find("H") !=
                 string::npos) { // current controlled voltage source
        valid = IsOkParamteres(words, 6);
        if (!valid) {
          cerr << "C.C.V.S. has invalid parameters " << endl;
          inputFile.close();
          return false;
        } else {
          FillDependantData(words, 0);
          numOfH++;
        }
      } else if (words[0].find("F") !=
                 string::npos) { // current controlled current source
        valid = IsOkParamteres(words, 6);
        if (!valid) {
          cerr << "C.C.C.S. has invalid parameters " << endl;
          inputFile.close();
          return false;
        } else {
          FillDependantData(words, 1);
          numOfF++;
        }
      } else if (words[0].find("E") !=
                 string::npos) { // voltage controlled voltage source
        valid = IsOkParamteres(words, 6);
        if (!valid) {
          cerr << "V.C.V.S. has invalid parameters " << endl;
          inputFile.close();
          return false;
        } else {
          FillDependantData(words, 2);
          numOfE++;
        }
      } else if (words[0].find("G") !=
                 string::npos) { // voltage controlled current source
        valid = IsOkParamteres(words, 6);
        if (!valid) {
          cerr << "V.C.C.S. has invalid parameters " << endl;
          inputFile.close();
          return false;
        } else {
          FillDependantData(words, 3);
          numOfG++;
        }
      }
    }
  }
  inputFile.close();
  return true;
}

Component SearchByName(const string &name) {
  Component requested;
  for (int i = 0; i < components.size(); i++) {
    if (components[i].name == name) {
      requested = components[i];
      break;
    }
  }
  return requested;
}

void SearchByNodes(const int &node1, const int &node2,
                   vector<Component *> &parrallelComp) {
  for (int i = 0; i < components.size(); i++) {
    if (components[i].node1 == node1 && components[i].node2 == node2) {
      parrallelComp.push_back(&components[i]);
    } else if (components[i].node1 == node2 && components[i].node2 == node1) {
      parrallelComp.push_back(&components[i]);
    }
  }
}

void AddZeroVs(const int &node1, const int &node2, const int &num) {
  Component zeroVs;
  zeroVs.name = "Vsz" + to_string(num);
  names.push_back(zeroVs.name);
  zeroVs.node1 = node1;
  zeroVs.node2 = node2;
  zeroVs.value = 0;
  components.push_back(zeroVs);
}

int CalcNumOfNodes() {
  int max = 0;
  for (int i = 0; i < components.size(); ++i) {
    if (components[i].node1 > max)
      max = components[i].node1;
    if (components[i].node2 > max)
      max = components[i].node2;
  }
  return max;
}

void AddCompInSeries(const int &negNode, const int &posNode, const int &nx,
                     const int &i) {
  vector<Component *> parallelComp;
  SearchByNodes(negNode, posNode, parallelComp);
  for (int k = 0; k < parallelComp.size(); k++) {
    bool foundH = parallelComp[k]->name.find("H") != string::npos;
    bool foundF = parallelComp[k]->name.find("F") != string::npos;
    bool foundE = parallelComp[k]->name.find("E") != string::npos;
    bool foundG = parallelComp[k]->name.find("G") != string::npos;
    if (!(foundH || foundF || foundE || foundG)) {
      if (parallelComp[k]->node1 == posNode) {
        parallelComp[k]->node1 = nx;
      } else if (parallelComp[k]->node2 == posNode) {
        parallelComp[k]->node2 = nx;
      }
    }
  }
  AddZeroVs(nx, posNode, i);
  for (auto &p : parallelComp) {
    p = nullptr;
    delete p;
  }
  parallelComp.clear();
}

/*zero voltage sources are numbered such like :
1 - first group of zero voltage sources are for CCVS and they numbered from 1 to
numOfH. 2 - second group of zero voltage sources are for CCCS and the numbered
from (numOfH + 1) to (numOfH + numOfF).
*/

void Modifications() {
  if (numOfH > 0) {
    numOfZeroVs = numOfH;
    for (int i = 1; i <= numOfZeroVs; i++) {
      numOfNodes = CalcNumOfNodes();
      int nx = numOfNodes + 1;
      string name = "H" + to_string(i);
      Component CCVS = SearchByName(name);
      AddCompInSeries(CCVS.negNode, CCVS.posNode, nx, i);
    }
  }
  if (numOfF > 0) {
    numOfZeroVs += numOfF;
    int countF = 1;
    for (int i = (numOfH + 1); i <= numOfZeroVs; i++) {
      numOfNodes = CalcNumOfNodes();
      int nx = numOfNodes + 1;
      string name = "F" + to_string(countF);
      Component CCCS = SearchByName(name);
      AddCompInSeries(CCCS.negNode, CCCS.posNode, nx, i);
      countF++;
    }
  }
}

void Nodes(const int &numOfNodes) {
  vector<Component> node;
  Component empty;
  empty.name = "empty";
  bool NotVS, NotCS;
  for (int i = 1; i <= numOfNodes; i++) {
    for (int k = 0; k < components.size(); k++) {
      NotVS = components[k].name.find("Vs") == string::npos;
      NotCS = components[k].name.find("Cs") == string::npos;
      if (NotVS && NotCS) {
        if (components[k].node1 == i || components[k].node2 == i) {
          node.push_back(components[k]);
        } else {
          node.push_back(empty);
        }
      }
    }
    nodesVector.push_back(node);
    node.clear();
  }
}

complex<double> CalcYPerNode(const int &row) {
  complex<double> TotalY(0, 0);
  for (int col = 0; col < nodesVector[row].size(); col++) {
    if (nodesVector[row][col].name != "empty") {
      TotalY += nodesVector[row][col].Y;
    }
  }
  return TotalY;
}

complex<double> CalcYCommon(const int &row1, const int &row2) {
  complex<double> CommonY(0, 0);
  for (int col1 = 0; col1 < nodesVector[row1].size(); col1++) {
    if (nodesVector[row1][col1].name != "empty") {
      for (int col2 = 0; col2 < nodesVector[row2].size(); col2++) {
        if (nodesVector[row2][col2].name != "empty") {
          if (nodesVector[row1][col1].name == nodesVector[row2][col2].name) {
            CommonY += nodesVector[row2][col2].Y;
            break;
          }
        }
      }
    }
  }
  return (CommonY *= -1);
}

void NodesPolarity(const int &numOfVS) {
  vector<complex<double>> polarity;
  for (int i = 0; i < components.size(); i++) {
    if (components[i].name.find("Vs") != string::npos) {
      polarity.push_back(-(components[i].node1));
      polarity.push_back(components[i].node2);
    } else if (components[i].name.find("H") != string::npos) {
      polarity.push_back(-(components[i].node1));
      polarity.push_back(components[i].node2);
    } else if (components[i].name.find("E") != string::npos) {
      polarity.push_back(-(components[i].node1));
      polarity.push_back(components[i].node2);
    }
    if (polarity.size() != 0) {
      nodePolarity.push_back(polarity);
    }
    polarity.clear();
  }
}

void SetEntriesZeros(Eigen::MatrixXcd &source) {
  Eigen::MatrixXcd zeros = Eigen::MatrixXcd::Zero(source.rows(), source.cols());
  source << zeros;
}

void MatrixG(const int &numOfNodes, const int &numOfG) {
  Eigen::MatrixXcd g(numOfNodes, numOfNodes);
  SetEntriesZeros(g);
  for (int row = 0; row < numOfNodes; row++) {
    for (int col = 0; col < numOfNodes; col++) {
      if (row == col) {
        g(row, col) = CalcYPerNode(row);
      } else {
        g(row, col) = CalcYCommon(row, col);
      }
    }
  }
  if (numOfG > 0) {
    for (int i = 1; i <= numOfG; i++) {
      string Gname = "G" + to_string(i);
      Component VCCS = SearchByName(Gname);
      if (VCCS.node1 != 0) {
        if (VCCS.negNode != 0) {
          g(VCCS.node1 - 1, VCCS.negNode - 1) -= VCCS.value;
        }
        if (VCCS.posNode != 0) {
          g(VCCS.node1 - 1, VCCS.posNode - 1) += VCCS.value;
        }
      }
      if (VCCS.node2 != 0) {
        if (VCCS.negNode != 0) {
          g(VCCS.node2 - 1, VCCS.negNode - 1) += VCCS.value;
        }
        if (VCCS.posNode != 0) {
          g(VCCS.node2 - 1, VCCS.posNode - 1) -= VCCS.value;
        }
      }
    }
  }
  G = g;
}

void MatrixB(const int &numOfNodes, const int &numOfVS, const int &numOfH,
             const int &numOfF, const int &numOfE, const int &numOfZeroVs) {
  int numOfCols = numOfVS + numOfH + numOfZeroVs + numOfE;
  Eigen::MatrixXcd b(numOfNodes, numOfCols);
  SetEntriesZeros(b);
  for (int col = 0; col < numOfCols; col++) {
    complex<double> node(1, 0);
    for (int row = 0; row < numOfNodes; row++) {
      if (node == abs(nodePolarity[col][0])) {
        b(row, col) = nodePolarity[col][0] / abs(nodePolarity[col][0]);
      } else if (node == abs(nodePolarity[col][1])) {
        b(row, col) = nodePolarity[col][1] / abs(nodePolarity[col][1]);
      } else {
        b(row, col) = 0;
      }
      node.real(node.real() + 1);
    }
  }
  if (numOfF > 0) {
    for (int i = 1; i <= numOfF; i++) {
      int VszNum = numOfH + i;
      string VszName = "Vsz" + to_string(VszNum);
      Component CCCS = SearchByName("F" + to_string(i));
      for (int k = 0; k < names.size(); k++) {
        if (names[k] == VszName) {
          if (CCCS.node1 != 0) {
            b(CCCS.node1 - 1, k) += CCCS.value;
          }
          if (CCCS.node2 != 0) {
            b(CCCS.node2 - 1, k) -= CCCS.value;
          }
          break;
        }
      }
    }
  }
  B = b;
}

void MatrixC(const int &numOfH, const int &numOfF, const int &numOfE) {
  C = B.transpose();
  if (numOfF > 0) {
    for (int i = 1; i <= numOfF; i++) {
      int VszNum = numOfH + i;
      string VszName = "Vsz" + to_string(VszNum);
      Component CCCS = SearchByName("F" + to_string(i));
      for (int k = 0; k < names.size(); k++) {
        if (names[k] == VszName) {
          if (CCCS.node1 != 0) {
            C(k, CCCS.node1 - 1) -= CCCS.value;
          }
          if (CCCS.node2 != 0) {
            C(k, CCCS.node2 - 1) += CCCS.value;
          }
          break;
        }
      }
    }
  }
  if (numOfE > 0) {
    for (int i = 1; i <= numOfE; i++) {
      string Ename = "E" + to_string(i);
      Component VCVS = SearchByName(Ename);
      for (int k = 0; k < names.size(); k++) {
        if (names[k] == Ename) {
          if (VCVS.negNode != 0) {
            C(k, VCVS.negNode - 1) += VCVS.value;
          }
          if (VCVS.posNode != 0) {
            C(k, VCVS.posNode - 1) -= VCVS.value;
          }
          break;
        }
      }
    }
  }
}

void MatrixD(const int &numOfVS, const int &numOfH, const int &numOfE,
             const int &numOfZeroVs) {
  int numOfRows = numOfVS + numOfH + numOfZeroVs + numOfE;
  int numOfCols = numOfVS + numOfH + numOfZeroVs + numOfE;
  D = Eigen::MatrixXcd::Zero(numOfRows, numOfCols);
  for (int row = 0; row < names.size(); row++) {
    if (names[row].find("H") != string::npos) {
      stringstream ss(names[row]);
      char c;
      int num;
      ss >> c >> num;
      for (int col = 0; col < names.size(); col++) {
        string requestedName = "Vsz" + to_string(num);
        if (names[col] == requestedName) {
          Component CCVS = SearchByName(names[row]);
          D(row, col) -= CCVS.value;
          break;
        }
      }
    }
  }
}

void MatrixA(const int &numOfNodes, const int &numOfVS, const int &numOfH,
             const int &numOfF, const int &numOfE, const int &numOfZeroVs) {
  int n = numOfNodes + numOfVS + numOfH + numOfZeroVs + numOfE;
  Eigen::MatrixXcd a(n, n);
  SetEntriesZeros(a);
  MatrixG(numOfNodes, numOfG);
  MatrixB(numOfNodes, numOfVS, numOfH, numOfF, numOfE, numOfZeroVs);
  MatrixC(numOfH, numOfF, numOfE);
  MatrixD(numOfVS, numOfH, numOfE, numOfZeroVs);
  a << G, B, C, D;
  A = a;
}

void MatrixI(const int &numOfNodes) {
  Eigen::MatrixXcd i(numOfNodes, 1);
  SetEntriesZeros(i);
  map<int, complex<double>> NodesCS;
  map<int, complex<double>>::iterator it;
  for (int i = 0; i < components.size(); i++) {
    if (components[i].name.find("Cs") != string::npos) {
      NodesCS[components[i].node1] += -(components[i].current);
      NodesCS[components[i].node2] += (components[i].current);
    }
  }
  for (it = NodesCS.begin(); it != NodesCS.end(); ++it) {
    if (it->first != 0) {
      i(it->first - 1, 0) = it->second;
    }
  }
  I = i;
}

void MatrixE(const int &numOfVS, const int &numOfH, const int &numOfE,
             const int &numOfZeroVs) {
  int numOfRows = numOfVS + numOfH + numOfZeroVs + numOfE;
  Eigen::MatrixXcd e(numOfRows, 1);
  SetEntriesZeros(e);
  int node = 0;
  for (int i = 0; i < components.size(); i++) {
    if (components[i].name.find("Vs") != string::npos) {
      for (int k = 0; k < names.size(); k++) {
        if (names[k] == components[i].name) {
          node = k;
          break;
        }
      }
      e(node, 0) = components[i].voltage;
      node++;
    }
  }
  E = e;
}

void MatrixZ(const int &numOfNodes, const int &numOfVS, const int &numOfH,
             const int &numOfE, const int &numOfZeroVs) {
  int n = numOfNodes + numOfVS + numOfH + numOfZeroVs + numOfE;
  Eigen::MatrixXcd z(n, 1);
  SetEntriesZeros(z);
  MatrixI(numOfNodes);
  MatrixE(numOfVS, numOfH, numOfE, numOfZeroVs);
  z << I, E;
  Z = z;
}

void MatrixX(const int &numOfNodes, const int &numOfVS, const int &numOfH,
             const int &numOfE, const int &numOfZeroVs) {
  int n = numOfNodes + numOfVS + numOfH + numOfZeroVs + numOfE;
  Eigen::MatrixXcd x(n, 1);
  SetEntriesZeros(x);
  x += A.jacobiSvd(ComputeThinU | ComputeThinV).solve(Z);
  for (int i = 0; i < n; i++) {
    x(i, 0).real(int(x(i, 0).real() * 1000) / 1000.0);
    x(i, 0).imag(int(x(i, 0).imag() * 1000) / 1000.0);
  }
  X = x;
}

// setting the voltage for all components except the VS.
void SetVoltages() {
  complex<double> v1, v2;
  for (int i = 0; i < components.size(); i++) {
    bool CCVS = components[i].name.find("H") != string::npos;
    bool VCVS = components[i].name.find("E") != string::npos;
    if (components[i].name.find("Vs") == string::npos) {
      if (components[i].node1 == 0) {
        v1 = 0;
        v2 = X(components[i].node2 - 1, 0);
      } else if (components[i].node2 == 0) {
        v1 = X(components[i].node1 - 1, 0);
        v2 = 0;
      } else {
        v1 = X(components[i].node1 - 1, 0);
        v2 = X(components[i].node2 - 1, 0);
      }
      components[i].voltage = v1 - v2;
      /*if (CCVS || VCVS) {
              complex<double> vnc, vpc;
              if (components[i].negNode != 0) {
                      vnc = X(components[i].negNode - 1, 0);
              }
              else {
                      vnc = (0,0);
              }
              if (components[i].posNode != 0) {
                      vpc = X(components[i].posNode - 1, 0);
              }
              else {
                      vpc = (0, 0);
              }
              double vncMag, vpcMag;
              vncMag = sqrt(pow(vnc.real(), 2) + pow(vnc.imag(), 2));
              vpcMag = sqrt(pow(vpc.real(), 2) + pow(vpc.imag(), 2));
              if (vpcMag <= vncMag) {
                      components[i].voltage = abs(components[i].voltage);
              }
              else {
                      components[i].voltage = -1 * abs(components[i].voltage);
              }
      }*/
      if (CCVS || VCVS)
        components[i].voltage = abs(components[i].voltage);
    }
  }
}

// setting the currents through the VS.
void setCurrents() {
  int pos = X.size() - 1;
  for (int i = components.size() - 1; i >= 0; i--) {
    bool VS = components[i].name.find("Vs") != string::npos;
    bool CCVS = components[i].name.find("H") != string::npos;
    bool VCVS = components[i].name.find("E") != string::npos;
    if (VS || CCVS || VCVS) {
      components[i].current = X(pos);
      pos--;
    }
  }
  if (numOfF > 0) {
    int Fnum = 1;
    string VszName;
    Component Vsz;
    for (int i = 0; i < components.size(); i++) {
      if (components[i].name.find("F") != string::npos) {
        VszName = "Vsz" + to_string(numOfH + Fnum);
        Vsz = SearchByName(VszName);
        components[i].current = Vsz.current * components[i].value;
        Fnum++;
      }
    }
  }
  if (numOfG > 0) {
    complex<double> v1 = 0;
    complex<double> v2 = 0;
    complex<double> vx = 0;
    for (int i = 0; i < components.size(); i++) {
      if (components[i].name.find("G") != string::npos) {
        if (components[i].node1 == 0) {
          v1 = 0;
          v2 = X(components[i].node2 - 1, 0);
        } else if (components[i].node2 == 0) {
          v1 = X(components[i].node1 - 1, 0);
          v2 = 0;
        } else {
          v1 = X(components[i].node1 - 1, 0);
          v2 = X(components[i].node2 - 1, 0);
        }
        vx = v1 - v2;
        components[i].current = components[i].value * vx;
      }
    }
  }
}

// calculate complex power for each element.first for R ,L and C.
void CalcPowerRLC() {
  bool R, L, C;
  complex<double> VeffSquared;
  for (int i = 0; i < components.size(); i++) {
    R = components[i].name.find("R") != string::npos;
    L = components[i].name.find("L") != string::npos;
    C = components[i].name.find("C") != string::npos;
    if (R || L || C) {
      VeffSquared = pow(components[i].voltage, 2);
      if (VeffSquared.imag() != 0) {
        VeffSquared /= 2;
      }
      if ((w == 0) && (L || C)) {
        components[i].power = 0;
      } else {
        // if (R) VeffSquared = abs(VeffSquared);
        components[i].power = VeffSquared / components[i].Zconj;
      }
    }
  }
}

// calculate the complex power for VS.
void CalcPowerVS() {
  complex<double> Veff;
  complex<double> IeffConj;
  for (int i = 0; i < components.size(); i++) {
    bool VS = components[i].name.find("Vs") != string::npos;
    bool CCVS = components[i].name.find("H") != string::npos;
    bool VCVS = components[i].name.find("E") != string::npos;
    if (VS || CCVS || VCVS) {
      Veff = components[i].voltage;
      if (Veff.imag() != 0) {
        Veff /= sqrt(2);
      }
      IeffConj = conj(components[i].current);
      if (IeffConj.imag() != 0) {
        IeffConj /= sqrt(2);
      }
      components[i].power = Veff * IeffConj;
    }
  }
}

// calculate the complex power for CS.
void CalcPowerCS() {
  complex<double> Veff;
  complex<double> IeffConj;
  for (int i = 0; i < components.size(); i++) {
    bool CS = components[i].name.find("Cs") != string::npos;
    bool CCCS = components[i].name.find("F") != string::npos;
    bool VCCS = components[i].name.find("G") != string::npos;
    if (CS || CCCS || VCCS) {
      Veff = components[i].voltage;
      if (Veff.imag() != 0) {
        Veff /= sqrt(2);
      }
      IeffConj = conj(components[i].current);
      if (IeffConj.imag() != 0) {
        IeffConj /= sqrt(2);
      }
      components[i].power = Veff * IeffConj;
    }
  }
}

bool DisplayX() {
  outputFile.open("output.txt", ios::trunc);
  if (!outputFile.is_open()) {
    cerr << "Error hapened while open the output file. " << endl;
    outputFile.close();
    return false;
  } else {
    outputFile << "==========unknowns vector================" << endl;
    outputFile << "[ ";
    for (int i = 0; i < X.size(); i++) {
      if (i != X.size() - 1) {
        outputFile << X(i, 0) << " , ";
      } else {
        outputFile << X(i, 0);
      }
    }
    outputFile << " ]" << endl;
    outputFile.close();
    return true;
  }
}

bool DisplayPower() {
  outputFile.open("output.txt", ios::app);
  if (!outputFile.is_open()) {
    cerr << "Error hapened while open the output file. " << endl;
    outputFile.close();
    return false;
  } else {
    outputFile << "\n==========power vector================" << endl;
    outputFile << "[ " << endl;
    for (int i = 0; i < components.size(); i++) {
      if (components[i].name.find("R") != string::npos) {
        double p = sqrt(pow(components[i].power.real(), 2) +
                        pow(components[i].power.imag(), 2));
        outputFile << components[i].name << " complex power is " << p << endl;
      } else {
        outputFile << components[i].name << " complex power is "
                   << components[i].power << endl;
      }
    }
    outputFile << "]" << endl;
    outputFile.close();
    return true;
  }
}

void Solver() {
  if (!ReadCircuitData()) {
    return;
  } else {
    Modifications();
    numOfNodes = CalcNumOfNodes();
    Nodes(numOfNodes);
    NodesPolarity(numOfVS);
    MatrixA(numOfNodes, numOfVS, numOfH, numOfF, numOfE, numOfZeroVs);
    MatrixZ(numOfNodes, numOfVS, numOfH, numOfE, numOfZeroVs);
    MatrixX(numOfNodes, numOfVS, numOfH, numOfE, numOfZeroVs);
    // for testing each matrix
    /*cout << A << "\n\n\n";
    cout << Z << "\n\n\n";
    cout << X;
    cout << "\n";*/
    SetVoltages();
    setCurrents();
    CalcPowerRLC();
    CalcPowerVS();
    CalcPowerCS();
    DisplayX();
    DisplayPower();
  }
}

int main() {
  Solver();
  cout << "Done." << endl;
  system("pause");
  return 0;
}
