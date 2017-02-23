//---------------------------------*-SWIG-*----------------------------------//
/*!
 * \file   parameterlist/Teuchos_XML.i
 * \author Seth R Johnson
 * \date   Tue Dec 06 18:01:54 2016
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
%{
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
%}

%inline %{
void load_from_xml(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const char* STRING, int SIZE)
{
    Teuchos::updateParametersFromXmlFile(std::string(STRING, SIZE),
                                         Teuchos::inOutArg(*plist));
}

void save_to_xml(const Teuchos::ParameterList& plist,
                 const char* STRING, int SIZE)
{
    Teuchos::writeParameterListToXmlFile(plist, std::string(STRING, SIZE));
}
%}


//---------------------------------------------------------------------------//
// end of parameterlist/Teuchos_XML.i
//---------------------------------------------------------------------------//
