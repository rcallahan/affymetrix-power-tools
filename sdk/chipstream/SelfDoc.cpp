////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#include "chipstream/SelfDoc.h"

bool SelfDoc::Opt::checkOption() const {
      bool minOk = false, maxOk = false;
      bool doThrow = Err::getThrowStatus();
      Err::setThrowStatus(true); 
      try {
        switch (type) {
        case SelfDoc::Opt::String :
            // no check for strings...
            minOk = true;
            maxOk = true;
            break;
        case SelfDoc::Opt::Integer :
            minOk = minVal == "NA" || (asInt() >= Convert::toInt(minVal.c_str()));
            maxOk = maxVal == "NA" || (asInt() <= Convert::toInt(maxVal.c_str()));
            break;
        case SelfDoc::Opt::Float :
            minOk = minVal == "NA" || (asFloat() >= Convert::toFloat(minVal.c_str()));
            maxOk = maxVal == "NA" || (asFloat() <= Convert::toFloat(maxVal.c_str()));
            break;
        case SelfDoc::Opt::Double :
            minOk = minVal == "NA" || (asDouble() >= Convert::toDouble(minVal.c_str()));
            maxOk = maxVal == "NA" || (asDouble() <= Convert::toDouble(maxVal.c_str()));
            break;
        case SelfDoc::Opt::Boolean :
            minOk = value == "true" || value == "false";
            maxOk = value == "true" || value == "false";
            break;
        default:
            Err::errAbort("SelfDoc::Opt::checkOption() - Don't recognize type: " + ToStr(type));
        }
      }
      catch(Except &e) {
        Err::setThrowStatus(doThrow);
        Err::errAbort("SelfDoc::Opt::checkOption() - Invalid value '"+name+"' specified for option '"+value+"'. Error: " + e.what());
      }
      catch(...) {
        Err::setThrowStatus(doThrow);
        Err::errAbort("SelfDoc::Opt::checkOption() - Invalid value '"+name+"' specified for option '"+value+"'.");
      }
      Err::setThrowStatus(doThrow);

      return minOk && maxOk;
    }

double SelfDoc::Opt::asDouble() const { return Convert::toDouble(value.c_str()); }
float SelfDoc::Opt::asFloat() const { return Convert::toFloat(value.c_str()); }
int SelfDoc::Opt::asInt() const { return Convert::toInt(value.c_str()); }
bool SelfDoc::Opt::asBool() const { return Convert::toBool(value.c_str()); }
std::string SelfDoc::Opt::asString() const { return value; }

std::string SelfDoc::getState() {
    std::vector<SelfDoc::Opt>::iterator vecIx;
    std::string description = getDocName();
    for(vecIx = m_Options.begin(); vecIx != m_Options.end(); ++vecIx) {
        description += ".";
        description += vecIx->name + "=" + vecIx->asString();
    }
    return description;
}

SelfDoc::Opt & SelfDoc::getDocOption(const std::string &name) { 
    std::map<std::string, int>::iterator i;
    i = m_OptionsIndex.find(name);
    if(i == m_OptionsIndex.end()) {
      Err::errAbort("SelfDoc::getDocOption() - Can't find option for name: " + ToStr(name));
    }
    return *(m_Options.begin() + i->second); 
  }

void SelfDoc::setOptValue(const std::string &name, const std::string &value) {
    SelfDoc::Opt &opt = getDocOption(name);
    opt.value = value;
    if(!opt.checkOption()) {
        Err::errAbort("SelfDoc::setOptValue() - '" + value + "' is not a valid option for '" + name + "'.");
    }
}

void SelfDoc::setOptValue(const std::string &name, const bool &value) {
    SelfDoc::Opt &opt = getDocOption(name);
    opt.value = value ? "true" : "false";
    if(!opt.checkOption()) {
        Err::errAbort("SelfDoc::setOptValue() - '" + ToStr(value ? "true" : "false") + 
                      "' is not a valid option for " + name);
    }
}

void SelfDoc::setDocOptions(const std::vector<SelfDoc::Opt> &param) { 
    std::map<std::string, int>::iterator offset;
    std::vector<SelfDoc::Opt>::const_iterator i;
    int count = 0;
    m_OptionsIndex.clear();
    m_Options.clear();
    for(i = param.begin(); i != param.end(); i++) {
        offset = m_OptionsIndex.find(i->name);
      if(offset != m_OptionsIndex.end()) {
          Err::errAbort("SelfDoc::setDocOptions() - Option name '" + i->name + 
                        "' has already been seen, names must be unique.");
      }
      if(!i->checkOption()) {
          Err::errAbort("SelfDoc::setDocOptions() - Option '" + i->name + 
                        "' doesn't pass its checks.");
      }
      m_OptionsIndex[i->name] = count++;
      m_Options.push_back(*i);
    }
}

SelfDoc SelfDoc::explainSelf() {
    SelfDoc doc;
    doc.setDocName("no name set");
    doc.setDocDescription("Currently undocumented.");
    // no parameters.
    return doc;
}

void SelfDoc::printExplanation(SelfDoc &doc, std::ostream &out) {
    std::string name = doc.getDocName();
    std::string description = doc.getDocDescription();
    std::vector<SelfDoc::Opt> opts = doc.getDocOptions();
    unsigned int maxLength = 0;
    std::vector<SelfDoc::Opt>::iterator iter;
    out << name << ":\n";
    Util::printStringWidth(out, description.c_str(), 0, name.size() + 2);
    out << std::endl;
    for(iter = opts.begin(); iter != opts.end(); iter++) {
        if(iter->name.size() > maxLength) 
            maxLength = iter->name.size();
    }
    maxLength += 3; // add characters for spaces and quotes
    
    out << std::endl << "Parameters: " << std::endl;    
    if(opts.empty()) {
        out << "   --- No Parameters ---   " << std::endl;
    }
    else {
        for(iter = opts.begin(); iter != opts.end(); iter++) {
            unsigned int currentLength = 0;
            std::string param = iter->name;
            std::string explanation = iter->descript + 
                " [default '" + iter->defaultVal+ "']";
            out << " '" << param << "' ";
            currentLength += param.size();
            while(currentLength < maxLength - 3) {
                out.put(' ');
                currentLength++;
            }
            Util::printStringWidth(out, explanation.c_str(), maxLength+1, currentLength);
            out << std::endl;
        }
    }
}
