#ifndef TYPESWITCHER_H_
#define TYPESWITCHER_H_

#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include "vcf.h"
#include "hts.h"

using namespace std;

class TypeSwitcher {
private:
	template<typename Tbcf, Tbcf (*convert)(const uint8_t*), Tbcf missing, Tbcf vector_end>
	void get_read(vector<int32_t>& v);

	//HX: enable reading string for fields like FT
	void get_read_string(string& s);

public:
	bcf_fmt_t* fmt;
	uint8_t* fmt_p;

	void (TypeSwitcher::*read)(vector<int32_t>&);
	//HX
	void (TypeSwitcher::*read_string)(string&);

	TypeSwitcher();
	virtual ~TypeSwitcher() noexcept;

	void init(bcf_fmt_t* fmt) throw (runtime_error);
};

#endif
