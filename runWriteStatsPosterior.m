%
function runWriteStatsPosterior(varargin)
  plotTypeList = ["Parameter","ECH4","SourceFraction","Loss"];

  % fixed period
  fixedPeriods = [
    1750 2015; 1850 2015; 1940 1950; 1986 2000; 2003 2012;
    1700 1750; 1750 1900; 1900 1950; 1950 1960; 1960 1970;
    1970 1980; 1980 1990; 1990 1995; 1995 2000; 2000 2005;
    2005 2010; 2010 2015; 1991 1998; 1999 2006; 2007 2014;
    1992 1998; 1999 2005; 1999 2005; 2006 2012; 1750 1850;
    1993 1999; 2000 2006; 2007 2013
  ];
  
  % one year period
  halfYearPeriods = (1970.5:1:2015.5)';
  halfYearPeriods = [halfYearPeriods, halfYearPeriods];
  
  % combine period data
  periodList = {[fixedPeriods;halfYearPeriods]};

  % choose specific period you want to calculate period mean difference under multi scenario mean
  periodListDiff = {
    [1991 1998; 1999 2006],[1999 2006; 2007 2014],[1992 1998; 1999 2005],...
    [1999 2005; 2006 2012],[1750 1850; 2003 2012],[1993 1999; 2000 2006],...
    [2000 2006; 2007 2013]};

  if ~isempty(varargin)
    workdirList = varargin{1};
  else
    workdirList = pwd; 
  end

  model = Model(workdirList(1),"Posterior",9);

  for plotType = plotTypeList  
    fieldList = getFieldList(plotType,model.nameList);
    data = Data(workdirList,model,plotType,fieldList);
    
    if length(string(workdirList)) == 1 %single result
      data.WriteStatics('periodList',periodList);
    else %multi scenario mean
      data.WriteStaticsMultiScenarioHistMean('calcType',"PeriodMean",'periodList',periodList);
      data.WriteStaticsMultiScenarioHistMean('calcType',"PeriodDiff",'periodList',periodListDiff);
    end
  end
 
end

%%
function fieldList = getFieldList(plotType,nameList)
  switch plotType
    case 'Atmosphere'
      fieldList = nameList.atmos;
    case 'Parameter'
      fieldList = nameList.para;
      if ~ismember("life",nameList.para)
        fieldList = [fieldList,"life","d13Cbio","d13Cff","dDbio","dDff"]; 
      end 
    case {'ECH4'} 
      fieldList = ["tot","bio","bb","ff","anth_ff","natr_ff",...
        "anth_bio","natr_bio","natr","anth"]; 
    case {'Loss'} 
      fieldList = ["totLoss","oh","ohLife"]; 
    case {'SourceFraction'} 
      fieldList = ["bio","ff","bb","anth_ff","geo"];
    otherwise
      error(strcat("No such plotType: ",plotType))
  end
end