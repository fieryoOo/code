#ifndef BINT_H
#define BINT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctype.h>

class bint {
public :
	// constructors
	bint() : isPositive(true) {
	}

	bint(long long val) {
		isPositive = val >= 0;
		if (!isPositive) val *= -1;
		//num.clear();
		if (val == 0) num = std::vector<int>(1, 0);
		while (val) {
			num.push_back(val % 10);
			val /= 10;
		}
	}

	bint(std::string str) {
		// remove spaces
		str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
		// remove leading zeros
		size_t digstart = isdigit(str[0]) ? 0 : 1;
		str.erase(digstart, str.find_first_not_of('0', digstart)-digstart);
		// find 'e/E' position in scientific notation
		size_t epos = str.find_first_of("eE");
		if( epos == std::string::npos ) epos = str.size();
		// start parsing
		//bool isPositive = true;
		//std::vector<unsigned short> absV;
		if( epos == 0 ) {
			num = {0};
		} else {
			// parse int part (before e/E)
			std::string str_p1 = getIntStr( str, 0, epos, isPositive );
			for ( const auto c : str_p1 )
				num.push_back( c - '0' );
			// parse scientific part (after e/E) and make correction
			if( epos+1 < str.size() ) {
				bool isPos2;
				std::string str_p2 = getIntStr( str, epos+1, str.size(), isPos2 );
				int order = std::stoi( str_p2 );
				if( isPos2 ) {
					num.insert( num.end(), order, 0 );
				} else if( order > num.size() ) {
					num = {0};
				} else {
					num.erase( num.end()-order, num.end() );
				}
			}
			std::reverse(num.begin(), num.end());
		}
	}

	friend const bint operator+(const bint& n1, const bint& n2){
		bint result, first(n1), second(n2);

		if(n1.isPositive != n2.isPositive){
			first.isPositive = !first.isPositive;
			result = first - second;

			result.isPositive = (n1.Abs() > n2.Abs() ? n1.isPositive : n2.isPositive);
			return result;
		}

		int carry = 0;
		for(size_t i=0; i < std::max(first.num.size(), second.num.size()); i++) {
			int calc = first.get(i) + second.get(i) + carry;
			result.num.push_back(calc % 10);
			carry = calc / 10;
		}
		if (carry)
			result.num.push_back(carry);

		result.isPositive = n1.isPositive;
		return result;
	}
	void operator+=(const bint& n2) { *this = *this + n2; }

	friend const bint operator-(const bint& n1, const bint& n2){
		bint result, first(n1), second(n2);
		result.isPositive = n1.isPositive;

		if(n1.isPositive != n2.isPositive){
			first.isPositive = !first.isPositive;
			result = first + second;
			result.isPositive = n1.isPositive;

			return result;
		}

		if(n1.Abs() < n2.Abs()){
			first = n2;
			second = n1;
			result.isPositive = !result.isPositive;
		}

		int carry = 0;
		for(size_t i=0; i < std::max(first.num.size(), second.num.size()); ++i) {
			int calc = first.get(i) - second.get(i) + carry;
			carry = 0;
			if(calc < 0){
				carry = -1;
				calc += 10;
			}
			result.num.push_back(calc);
		}
		if (carry)
			result.num.push_back(carry);


		for(size_t i=result.num.size()-1; !result.num[i] && i>0; i--){
			result.num.pop_back();
		}

		return result;
	}
	void operator-=(const bint& n2) { *this = *this - n2; }

	friend const bint operator*(const bint& n1, const bint& n2){
		/*bint result, greater, less;
		greater = std::max(n1.Abs(), n2.Abs());
		less = std::min(n1.Abs(), n2.Abs());

		for(int i=0; i<less; i++){
			result = result + greater;
		}
		result.isPositive = (n1.isPositive == n2.isPositive ? true : false);

		if(result.num.empty()) result = 0;

		return result;*/

		bint result;

		int carry = 0;
		for(size_t i=0; i<n2.num.size(); i++){
			carry = 0;
			for(size_t j=0; j<n1.num.size(); j++){
				int calc = n2.get(i) * n1.get(j) + carry;
				carry = calc / 10;
				if(i+j >= result.num.size()){
					result.num.push_back(calc % 10);
				}
				else{
					result.num[i+j] += calc % 10;
					if(result.num[i+j] >= 10){
						carry += result.num[i+j] / 10;
						result.num[i+j] = result.num[i+j] % 10;
					}
				}
			}
			if(carry)
				result.num.push_back(carry);
		}

		for(size_t i=result.num.size()-1; !result.num[i] && i>0; i--){
			result.num.pop_back();
		}

		result.isPositive = (n1.isPositive == n2.isPositive ? true : false);
		if(result.Abs() == 0) result.isPositive = true;
		return result;
	}
	void operator*=(const bint& n2) { *this = *this * n2; }

	friend const bint operator/(const bint& n1, const bint& n2){
		if(n2 == 0)
			throw std::runtime_error("ArithmeticException: divided by zero");
		bint result(0);
		bint dividend(n1.Abs()), divisor(n2.Abs());

		while(dividend >= divisor){
			dividend = dividend - divisor;
			result = result + 1;
		}

		result.isPositive = (n1.isPositive == n2.isPositive);
		return result;
	}
	void operator/=(const bint& n2) { *this = *this / n2; }

	friend const bint operator%(const bint& n1, const bint& n2){
		if(n2 == 0)
			throw std::runtime_error("ArithmeticException: mod zero");

		bint result, share;

		share = n1 / n2;
		result = n1 - share * n2;

		return result;
	}
	void operator%=(const bint& n2) { *this = *this % n2; }

	friend bool operator==(const bint& left, const bint& right) {
		return left.isPositive == right.isPositive && left.num == right.num;
	}

	friend bool operator>(const bint& left, const bint& right){
		if(left.isPositive != right.isPositive){
			return left.isPositive ? true : false;	
		}

		if(left.num.size() > right.num.size()){
			return left.isPositive ? true : false;
		}
		else if(left.num.size() < right.num.size()){
			return left.isPositive ? false : true;
		}

		for(int i=left.num.size()-1; i>=0; --i){
			if(left.num[i] > right.num[i]) return left.isPositive ? true : false;
			else if(left.num[i] < right.num[i]) return left.isPositive ? false : true;
		}

		return false;
	}

	friend bool operator>=(const bint& left, const bint& right){
		return (left > right) || (left == right);
	}

	friend bool operator<(const bint& left, const bint& right){
		return !(left >= right);
	}

	friend bool operator<=(const bint& left, const bint& right){
		return (left < right || left == right);
	}

	friend bool operator!=(const bint& left, const bint& right) {
		return !(left.isPositive == right.isPositive && left.num == right.num);
	}

	friend const bint operator<<(const bint& srcNum, const bint& shiftCount){
		if(shiftCount < 0) 
			throw "shiftCount is not positive" ;
		bint result(srcNum);

		for(int i=0; i<shiftCount; i++){
			result = result * 2;
		}

		return result;
	}

	friend const bint operator>>(const bint& srcNum, const bint& shiftCount){
		if(shiftCount < 0) 
			throw "shiftCount is not positive";
		bint result(srcNum);

		for(int i=0; i<shiftCount; i++){
			if(srcNum < 0 && result.num[0] & 1) result = result - 1;
			result = result / 2;
		}

		return result;
	}

	friend bool operator&&(const bint& n1, const bint& n2){
		if(n1 == 0 || n2 == 0) return false;
		return true;
	}

	friend bool operator||(const bint& n1, const bint& n2){
		if(n1 == 0 && n2 == 0) return false;
		return true;
	}

	friend const bint operator&(const bint& n1, const bint& n2){
		bint result = n1.innerBitwise(n2, '&');
		return result;
	}

	friend const bint operator|(const bint& n1, const bint& n2){
		bint result = n1.innerBitwise(n2, '|');
		return result;
	}

	friend const bint operator^(const bint& n1, const bint& n2){
		bint result = n1.innerBitwise(n2, '^');
		return result;
	}

	const bint Abs() const{
		bint result(*this);
		result.isPositive = true;
		return result;
	}

	friend std::ostream& operator<<(std::ostream& os, const bint& val) { 
		if(!val.isPositive) os << "-";
		for(int i=val.num.size()-1; i>=0; --i){
			os << val.num[i];
		}
		return os;
	}

	void print() { std::cout<<this; }


private :
	std::vector<int> num;
	bool isPositive;

	inline int get(int digit) const {
		return num.size() <= digit ? 0 : num[digit];
	}

	inline std::string getIntStr( const std::string& str, const size_t is, const size_t ie, bool& isPos ) {
		// zero/negative length
		isPos = true;
		if( is >= ie ) return "0";
		size_t istart = isdigit(str[is]) ? is : is+1;
		// check sign character
		if( istart==is || str[is]=='+' ) {
			isPos = true;
		} else if( str[is] == '-' ) {
			isPos = false;
		} else {
			throw std::runtime_error("invalid int string");
		}
		// empty digit string
		if( istart == ie ) return "0";
		// search for invalid digits
		bool allzeros = true;
		for(int i=istart; i<ie; i++) {
			if( ! isdigit(str[i]) )
				throw std::runtime_error("invalid int string");
			if( str[i] != '0' ) allzeros = false;
		}
		// return 0 if all zeros
		if( allzeros ) return "0";
		return str.substr(istart, ie-istart);
	}

	inline const bint convertToBinary(){
		bint result;

		while(1){
			int calc = this->get(0) % 2;
			result.num.push_back(calc);
			*this = *this >> 1;
			if(*this == 0 || *this == -1) break;
		}
		
		return result;
	}

	inline const bint convertToDecimal(){
		bint result(0), digit = 1;
		
		for(size_t i=0; i<num.size(); i++){
			result = result + get(i) * digit;
			digit = digit * 2;
		}

		if(!isPositive){
			result = result - digit;
		}

		return result;
	}

	const bint innerBitwise(const bint& other, char op) const{
		bint result, n1_binary(*this), n2_binary(other);
		
		n1_binary = n1_binary.convertToBinary();
		n2_binary = n2_binary.convertToBinary();
		
		for(size_t i=0; i<std::min(n1_binary.num.size(), n2_binary.num.size()); i++){
			if(op == '&')
				result.num.push_back(n1_binary.get(i) & n2_binary.get(i));
			else if(op == '|')
				result.num.push_back(n1_binary.get(i) | n2_binary.get(i));
			else if(op == '^')
				result.num.push_back(n1_binary.get(i) ^ n2_binary.get(i));
			else
				throw std::runtime_error("Operator not define");
		}
		
		if(n1_binary.num.size() < n2_binary.num.size()){
			if(op == '&' && !this->isPositive){
				for(size_t i=n1_binary.num.size(); i<n2_binary.num.size(); i++){
					result.num.push_back(n2_binary.get(i));
				}
			}
			else if(op == '|' && this->isPositive){
				for(size_t i=n1_binary.num.size(); i<n2_binary.num.size(); i++){
					result.num.push_back(n2_binary.get(i));
				}
			}
			else if(op == '^'){
				for(size_t i=n1_binary.num.size(); i<n2_binary.num.size(); i++){
					int bit = n2_binary.get(i);
					if(!this->isPositive) bit = !bit;
					result.num.push_back(bit);
				}
			}
		}
		else if(n1_binary.num.size() > n2_binary.num.size()){
			if(op == '&' && !other.isPositive){
				for(size_t i=n2_binary.num.size(); i<n1_binary.num.size(); i++){
					result.num.push_back(n1_binary.get(i));
				}
			}
			else if(op == '|' && other.isPositive){
				for(size_t i=n2_binary.num.size(); i<n1_binary.num.size(); i++){
					result.num.push_back(n1_binary.get(i));
				}
			}
			else if(op == '^'){
				for(size_t i=n2_binary.num.size(); i<n1_binary.num.size(); i++){
					int bit = n1_binary.get(i);
					if(!other.isPositive) bit = !bit;
					result.num.push_back(bit);
				}
			}
		}
		
		if(op == '&'){
			if(!this->isPositive && !other.isPositive)
				result.isPositive = false;
		}
		else if(op == '|'){
			if(!this->isPositive || !other.isPositive)
				result.isPositive = false;
		}
		else if(op == '^'){
			if(this->isPositive != other.isPositive)
				result.isPositive = false;
		}
		
		result = result.convertToDecimal();
		return result;
	}

};

#endif
