import re
import os

objects = ['Photon', 'Electron', 'Muon', 'Jet', 'Vertex']
susyObjects = {'Photon': 'Photon', 'Electron': 'Electron', 'Muon': 'Muon', 'Jet': 'PFJet', 'Vertex': 'Vertex'}

objectVars = file('ObjectVars.h')

classPat = re.compile('^[ ]*class[ ]+([a-zA-Z0-9]+)Vars[ ]*{')
cTorPat = re.compile('^[ ]*[a-zA-Z0-9]+Vars\([^,]+(,[ ]+Event.*|)\);')
varPat = re.compile('^[ ]*((?:unsigned[ ]|)(?:bool|char|short|int|unsigned|long|float|double))[ ]+([a-zA-Z_][a-zA-Z0-9_]*);')

useEvent = dict()
varList = dict()

obj = ''
for line in objectVars:
    if '};' in line:
        obj = ''
        
    if obj:
        cTorMatch = cTorPat.match(line)
        if cTorMatch:
            useEvent[obj] = len(cTorMatch.group(1)) != 0
            
        varMatch = varPat.match(line)
        if varMatch:
            varList[obj].append((varMatch.group(1), varMatch.group(2)))
            
    lineMatch = classPat.match(line)
    if lineMatch and lineMatch.group(1) in objects:
        obj = lineMatch.group(1)
        varList[obj] = []

objectVars.close()

# GENERATE HEADER

headerContent = '''/* Auto-generated header file */
#ifndef ObjectTree_h
#define ObjectTree_h

#include "ObjectVars.h"

#include "TTree.h"
#include "TString.h"

namespace susy {

  unsigned const NMAX(512);
'''
 
for obj in objects:
    headerContent += '''
  class ''' + obj + '''VarsArray {
  public:
    ''' + obj + '''VarsArray() : size(0) {}
    ~''' + obj + '''VarsArray() {}
    void setBranches(TTree&);
    void setAddress(TTree&);
    void push_back(''' + obj + '''Vars const&);
    void clear() { size = 0; }
    ''' + obj + '''Vars at(unsigned) const;

    unsigned size;
'''

    for (type, name) in varList[obj]:
        headerContent += '''
    ''' + type + ' ' + name + '[NMAX];'

    headerContent += '''
  };
'''

headerContent += '''
  class ObjectTree {
  public:    
    ObjectTree();
    ~ObjectTree();

    void setOutput(TString const&,'''

for i in range(len(objects)):
    headerContent += ' bool = true'
    if i != len(objects) - 1:
        headerContent += ','
    else:
        headerContent += ');'

headerContent += '''
    void setOutput(TTree&,'''

for i in range(len(objects)):
    headerContent += ' bool = true'
    if i != len(objects) - 1:
        headerContent += ','
    else:
        headerContent += ');'

headerContent += '''
    static void setBranchStatus(TTree&,'''

for i in range(len(objects)):
    headerContent += ' bool = true'
    if i != len(objects) - 1:
        headerContent += ','
    else:
        headerContent += ');'

headerContent += '''
    void initEvent(Event const&);
    void fill() { output_->Fill(); }'''

for obj in objects:
    lowerName = obj.lower()
    headerContent += '''
    void save(''' + obj + 'Vars const& _vars) { ' + lowerName + 'Array_.push_back(_vars); }'

for obj in objects:
    lowerName = obj.lower()
    headerContent += '''
    unsigned get''' + obj + 'Size() const { return ' + lowerName + 'Array_.size; }'

for obj in objects:
    lowerName = obj.lower()
    headerContent += '''
    ''' + obj + 'VarsArray const& get' + obj + 'Array() const { return ' + lowerName + 'Array_; }'

headerContent += '''
  private:
    void setBranches_('''

for i in range(len(objects)):
    headerContent += 'bool'
    if i != len(objects) - 1:
        headerContent += ', '
    else:
        headerContent += ');'

for obj in objects:
    lowerName = obj.lower()
    headerContent += '''
    ''' + obj + '''VarsArray ''' + lowerName + '''Array_;'''

headerContent += '''
    unsigned runNumber_;
    unsigned lumiNumber_;
    unsigned eventNumber_;

    TTree* output_;
    bool ownOutput_;
  };

}

#endif
'''

headerFile = file('ObjectTree.h', 'w')
headerFile.write(headerContent)
headerFile.close()

# GENERATE SRC

cTors = dict()
setBranches = dict()
setAddress = dict()
pushBack = dict()
at = dict()

for obj in objects:
    lowerName = obj.lower()

    cTorText = '''
  ''' + obj + 'Vars::' + obj + '''Vars() :'''

    initList = ''
    for (type, name) in varList[obj]:
        initList += '''
    ''' + name + '('
        if type == 'float' or type == 'double':
            initList += '0.'
        elif type == 'bool':
            initList += 'false'
        else:
            initList += '0'

        initList += '),'

    initList = initList.rstrip(',')
    cTorText += initList
    cTorText += '''
  {
  }
'''
    cTors[obj] = cTorText

    setBranchText = '''
  void
  ''' + obj + '''VarsArray::setBranches(TTree& _tree)
  {
    _tree.Branch("''' + lowerName + '.size", &size, "' + lowerName + '.size/i");'

    for (type, name) in varList[obj]:
        branch = '''
    _tree.Branch("''' + lowerName + '.' + name + '", ' + name + ', "' + name + '[' + lowerName + '.size]/'
        if type == 'char':
            branch += 'B'
        elif type == 'unsigned char':
            branch += 'b'
        elif type == 'short':
            branch += 'S'
        elif type == 'unsigned short':
            branch += 's'
        elif type == 'int':
            branch += 'I'
        elif type == 'unsigned' or type == 'unsigned int':
            branch += 'i'
        elif type == 'long':
            branch += 'L'
        elif type == 'unsigned long':
            branch += 'l'
        elif type == 'float':
            branch += 'F'
        elif type == 'double':
            branch += 'D'
        elif type == 'bool':
            branch += 'O'

        branch += '");'
        setBranchText += branch

    setBranchText += '''
  }
'''
    setBranches[obj] = setBranchText

    setAddressText = '''
  void
  ''' + obj + '''VarsArray::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    _tree.SetBranchAddress("''' + lowerName + '.size", &size);'

    for (type, name) in varList[obj]:
        bName = lowerName + '.' + name
        setAddressText += '''
    if(_tree.GetBranch("''' + bName + '")) _tree.SetBranchAddress("' + bName + '", ' + name + ''');
    else notFound.push_back("''' + bName + '");'

    setAddressText += '''
    
    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
  }
'''
    setAddress[obj] = setAddressText

    pushBackText = '''
  void
  ''' + obj + 'VarsArray::push_back(' + obj + '''Vars const& _vars)
  {
    if(size == NMAX - 1)
      throw Exception(Exception::kEventAnomaly, "Too many ''' + obj + '''s");
'''
    
    for (type, name) in varList[obj]:
        pushBackText += '''
    ''' + name + '[size] = _vars.' + name + ';'

    pushBackText += '''
    ++size;
  }
'''
    pushBack[obj] = pushBackText

    atText = '''
  ''' + obj + '''Vars
  ''' + obj + '''VarsArray::at(unsigned _pos) const
  {
    if(_pos >= size)
      throw std::range_error("''' + obj + '''Vars out-of-bounds");
      
    ''' + obj + '''Vars vars;
'''

    for (type, name) in varList[obj]:
        atText += '''
    vars.''' + name + ' = ' + name + '[_pos];'

    atText += '''
    return vars;
  }
'''

    at[obj] = atText

preamble = '#include "ObjectVars.h"\n'

try:
    originalSrc = file('ObjectVars.cc', 'r')
    userDef = ''

    copy = False
    namespace = False
    for line in originalSrc:
        if 'namespace susy' in line:
            namespace = True
            
        if not namespace and 'ObjectVars.h' not in line and not re.match('^[ ]*/\*.*\*/[ ]*$', line):
            preamble += line

        if '/* START USER-DEFINED IMPLEMENTATION (DO NOT MODIFY THIS LINE) */' in line:
            copy = True

        if copy:
            userDef += line

        if '/* END USER-DEFINED IMPLEMENTATION (DO NOT MODIFY THIS LINE) */' in line:
            copy = False

    originalSrc.close()
    
except:
    userDef = '\n/* START USER-DEFINED IMPLEMENTATION (DO NOT MODIFY THIS LINE) */\n'

    for obj in objects:
        userDef += '''
  void
  ''' + obj + '''Vars::set(''' + susyObjects[obj] + ' const&'
        if useEvent[obj]:
            userDef += ', Event const&'

        userDef += ''')
  {
  }

  /*static*/
  ''' + obj + '''Vars::setBranchStatus(TTree&)
  {
  }
'''

    userDef += '/* END USER-DEFINED IMPLEMENTATION (DO NOT MODIFY THIS LINE) */\n'

# ObjectTree.cc

objTreeContent = '''/* Auto-generated source file */
#include "ObjectTree.h"
#include "Utilities.h"

#include "TFile.h"

#include <iostream>
#ifndef STANDALONE
#include "SusyEvent.h"
#endif

namespace susy {
'''

for obj in objects:
    objTreeContent += setBranches[obj]
    objTreeContent += setAddress[obj]
    objTreeContent += pushBack[obj]
    objTreeContent += at[obj]

objTreeContent += '''

  ObjectTree::ObjectTree() :'''

for obj in objects:
    lowerName = obj.lower()
    objTreeContent += '''
    ''' + lowerName + '''Array_(),'''

objTreeContent += '''
    runNumber_(0),
    lumiNumber_(0),
    eventNumber_(0),
    output_(0),
    ownOutput_(false)
  {
  }

  ObjectTree::~ObjectTree()
  {
    if(ownOutput_ && output_){
      TFile* outFile(output_->GetCurrentFile());
      outFile->cd();
      output_->Write();
      delete outFile;
    }
  }

  void
  ObjectTree::setOutput(TString const& _fileName'''

for obj in objects:
    objTreeContent += ', bool _set' + obj + '/* = true*/'

objTreeContent += ''')
  {
    ownOutput_ = true;

    TFile::Open(_fileName, "recreate");
    output_ = new TTree("objectVars", "Object ID variables");

    setBranches_('''

for obj in objects:
    objTreeContent += '_set' + obj + ', '

objTreeContent = objTreeContent.rstrip(', ')
objTreeContent += ''');
  }

  void
  ObjectTree::setOutput(TTree& _tree'''

for obj in objects:
    objTreeContent += ', bool _set' + obj + '/* = true*/'

objTreeContent += ''')
  {
    output_ = &_tree;

    setBranches_('''

for obj in objects:
    objTreeContent += '_set' + obj + ', '

objTreeContent = objTreeContent.rstrip(', ')
objTreeContent += ''');
  }

  /*static*/
  void
  ObjectTree::setBranchStatus(TTree& _input'''

for obj in objects:
    objTreeContent += ', bool _set' + obj + '/* = true*/'

objTreeContent += ''')
  {
    _input.SetBranchStatus("runNumber", 1);
    _input.SetBranchStatus("luminosityBlockNumber", 1);
    _input.SetBranchStatus("eventNumber", 1);
'''

for obj in objects:
    objTreeContent += '''
    if(_set''' + obj + ') ' + obj + 'Vars::setBranchStatus(_input);'

objTreeContent += '''
  }

#ifdef STANDALONE
  void
  ObjectTree::initEvent(Event const&)
  {
    runNumber_ = 0;
    lumiNumber_ = 0;
    eventNumber_ = 0;
#else
  void
  ObjectTree::initEvent(Event const& _event)
  {
    runNumber_ = _event.runNumber;
    lumiNumber_ = _event.luminosityBlockNumber;
    eventNumber_ = _event.eventNumber;
#endif'''

for obj in objects:
    objTreeContent += '''
    ''' + obj.lower() + 'Array_.clear();'

objTreeContent += '''
  }

  void
  ObjectTree::setBranches_('''

for obj in objects:
    objTreeContent += 'bool _set' + obj + ', '

objTreeContent = objTreeContent.rstrip(', ') + ')'
objTreeContent += '''
  {
    output_->Branch("runNumber", &runNumber_, "runNumber/i");
    output_->Branch("lumiNumber", &lumiNumber_, "lumiNumber/i");
    output_->Branch("eventNumber", &eventNumber_, "eventNumber/i");
'''

for obj in objects:
    objTreeContent += '''
    if(_set''' + obj + ') ' + obj.lower() + 'Array_.setBranches(*output_);'

objTreeContent += '''
  }
'''

objTreeContent += '}\n'

objTreeFile = file('ObjectTree.cc', 'w')
objTreeFile.write(objTreeContent)
objTreeFile.close()

# ObjectVars.cc

objVarsContent = '''/* Partially auto-generated source file - edit where indicated */
/* Add necessary inclusions below */
''' + preamble + '''namespace susy {
'''

for obj in objects:
    objVarsContent += cTors[obj]

objVarsContent += '\n'
objVarsContent += userDef
objVarsContent += '''
}
'''

objVarsFile = file('ObjectVars.cc', 'w')
objVarsFile.write(objVarsContent)
objVarsFile.close()
