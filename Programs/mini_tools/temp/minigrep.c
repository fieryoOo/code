  /*
   * minigrep.cpp for (self)
   *
   * Made by tuxedo
   * Login <pierre@tuxedo.fr>
   *
   * Started on Wed 16 Mar 2007 15:18:39 tuxedo
   * Last update Fri 16 Mar 2007 14:22:37 tuxedo
   */
        
 #include <algorithm>
 #include <iostream>
 #include <string>
 #include <vector>
        
 using std::string;
 using std::vector;
      
 static string const WW_OPT = "-w";
 static string const AM_OPT = "-a";
        
 void displayLine(string const& s) {
     std::cout << s << std::endl;
 }
        
 bool lineIs(string const& line, string trigger) {
     // We must ignore ^ and $ here!
     if (!trigger.empty()) {
         if (^ == trigger[0])
             trigger.erase(trigger.begin());
         if (!trigger.empty() && $ == trigger[trigger.size() - 1])
             trigger.erase(trigger.end() - 1);
     }
     return line == trigger;
 } // lineIs
        
 bool lineContains(string const& line, string trigger) {
     if ("" == trigger)
         return true;
     if (^ == trigger[0]) {
         trigger.erase(trigger.begin());
         if (!trigger.empty() && $ == trigger[trigger.size() - 1])
             return lineIs(line, trigger);
         return 0 == line.find(trigger);
     }
     if ($ == trigger[trigger.size() - 1]) {
         trigger.erase(trigger.end() - 1);
         string::size_type const pos = line.rfind(trigger);
         return string::npos != pos && trigger.size() + pos == line.size();
     }
     return string::npos != line.find(trigger);
 } // lineContains
        
 void grepAdjacent(vector<string> const& lines, vector<string> const& triggers,
     bool matchWholeLines) {
         vector<string>::const_iterator it = lines.begin();
         while (lines.end() != (it = std::search(it, lines.end(), triggers.begin(), triggers.end(), matchWholeLines ? lineIs : lineContains))) {
       
       
             std::for_each(it, it + triggers.size(), displayLine);
             it += triggers.size();
         }
 } // grepAdjacent
        
 void grepDisjoint(vector<string> const& lines, vector<string> const& triggers, bool matchWholeLines) {
 
     vector<string>::const_iterator it = lines.begin();
     while (lines.end() != (it = std::find_first_of(it, lines.end(), triggers.begin(), triggers.end(), matchWholeLines ? lineIs : lineContains))) {
 
       
         displayLine(*it);
         ++it;
     }
 } // grepDisjoint
        
 int main(int argc, char const * argv[]) {
     // There is no simple way to force arguments at beginning of command line.
     // Therefore, the ugly, tortuous..
     bool matchWholeLines = 1 < argc && argv[1] == WW_OPT;
     bool useAdjacentMatches = matchWholeLines ? 2 < argc && argv[2] == AM_OPT : 1 < argc && argv[1] == AM_OPT;
 
       
     if (useAdjacentMatches && !matchWholeLines)
         matchWholeLines = 2 < argc && argv[2] == WW_OPT;
     // Now fetch search strings
     int const triggerStart = (matchWholeLines ? 1 : 0) + (useAdjacentMatches ? 1 : 0) + 1;
 
     vector<string> triggers;
     for (int index = triggerStart; index < argc; ++index)
         triggers.push_back(argv[index]);
     // Fetch standard input
     vector<string> lines;
     string line;
     while (getline(std::cin, line))
         lines.push_back(line);
     // Display results
     if (triggers.empty())
         std::for_each(lines.begin(), lines.end(), displayLine);
     else if (useAdjacentMatches)
        grepAdjacent(lines, triggers, matchWholeLines);
    else
        grepDisjoint(lines, triggers, matchWholeLines);
    return 0;
}
