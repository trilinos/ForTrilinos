/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
%}

%inline %{
void load_from_xml(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const std::string& xml_path)
{
    Teuchos::updateParametersFromXmlFile(xml_path, Teuchos::inOutArg(*plist));
}

void save_to_xml(const Teuchos::ParameterList& plist,
                 const std::string& xml_path)
{
    Teuchos::writeParameterListToXmlFile(plist, xml_path);
}
%}
