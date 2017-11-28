/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
%}

%inline %{
void load_from_xml(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   std::pair<const char*, size_t> xml)
{
    Teuchos::updateParametersFromXmlFile(std::string(xml.first, xml.second),
                                         Teuchos::inOutArg(*plist));
}

void save_to_xml(const Teuchos::ParameterList& plist,
                 std::pair<const char*, size_t> xml)
{
    Teuchos::writeParameterListToXmlFile(plist, std::string(xml.first, xml.second));
}
%}
