#ifndef FILENAME_H
#define FILENAME_H

#include "StackTrace.h"
#include <unistd.h>
#include <string>
#include <stdexcept>

#define FuncName __FUNCTION__

class FileName : public std::string {

   class ErrorBase : public std::runtime_error {
   public:
      ErrorBase(const std::string message)
			: runtime_error(message) {
			PrintStacktrace();
      }
   };
   class EmptyName : public ErrorBase {
   public:
      EmptyName(const std::string funcname, const std::string info = "")
         : ErrorBase("Error("+funcname+"): Empty file name ("+info+").") {}
   };
   class BadFile : public ErrorBase {
   public:
      BadFile(const std::string funcname, const std::string info = "")
			: ErrorBase("Error("+funcname+"): Cannot access file ("+info+").") {}
   };

public:
   FileName() : std::string() {}
   FileName(const char* fnamein) : std::string(fnamein) {}
   FileName(const std::string fnamein) : std::string(fnamein) {}

   void CheckEmpty() const {
      if( this->empty() ) throw EmptyName(FuncName);
   }
   void CheckAccess() const {
      if( this->empty() ) throw EmptyName(FuncName);
      if( access(this->c_str(), R_OK) == -1 ) throw BadFile(FuncName, *this);
   }

};

#endif
