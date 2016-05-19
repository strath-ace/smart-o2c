clear all; close all; clc
addpath(genpath(strcat(pwd,'/AsteroidMainBelt')));
addpath(genpath(strcat(fileparts(fileparts(pwd)),'/Optimisation/AIDMAP')));
addpath(strcat(fileparts(fileparts(pwd)),'/Optimisation'));

% This is the main file for the Atira problem
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Physarum Options         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.LinearDilationCoefficient = 5e-3;                       %Linear dilation coefficient 'm'
options.EvaporationCoefficient = 1e-4;                          %Evaporation coefficient 'rho'
options.GrowthFactorVal = 5e-3;                                 %Growth factor 'GF'
options.NumberOfAgents = 10;                                    %Number of virtual agents 'N_agents'
options.RamificationProbability = 0.7;                          %Probability of ramification 'p_ram'
options.RamificationWeight = 1;                                 %Weight on ramification 'lambda'
options.MaximumRadiusRatio = 20;                                %Maximum ratio between the link's radius & the starting radius
options.MinimumRadiusRatio = 1e-3;                              %Maximum ratio between the link's radius & the starting radius
options.StartingRadius = 1;                                     %The starting radius of the veins
options.RamificationAmount = 3;                                 %The number of nodes initially generated for the ramification
options.Generations = 40;                                       %The number of generations
options.Viscosity = 1;                                          %The viscocity of the "fluid" 
options.MinCommonNodesThres = 5;                                %The minimum number of nodes two decision sequences should have in common for a restart to occur
options.IfZeroLength = 1e-15;                                   %Value assigned to the length if it's zero (to prevent flux = inf)
options.MaxChildFindAttempts = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Problem-Specific Options     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.Targets = {'Ceres','Pallas','Juno','Vesta','Astraea','Hebe','Iris','Flora','Metis','Hygiea','Parthenope','Victoria','Egeria','Irene','Eunomia','Psyche','Thetis','Melpomene','Fortuna','Massalia','Lutetia','Kalliope','Thalia','Themis','Phocaea','Proserpina','Euterpe','Bellona','Amphitrite','Urania','Euphrosyne','Pomona','Circe','Leukothea','Atalante','Fides','Leda','Laetitia','Harmonia','Daphne','Isis','Ariadne','Nysa','Eugenia','Hestia','Aglaja','Doris','Pales','Virginia','Nemausa','Europa','Kalypso','Alexandra','Pandora','Melete','Mnemosyne','Concordia','Elpis','Echo','Danae','Erato','Ausonia','Maja','Asia','Leto','Hesperia','Panopaea','Niobe','Feronia','Galatea','Eurydike','Frigga','Diana','Eurynome','Sappho','Terpsichore','Alkmene','Beatrix','Klio','Io','Semele','Thisbe','Julia','Antiope','Aegina','Undina','Minerva','Aurora','Arethusa','Aegle','Klotho','Ianthe','Dike','Hekate','Helena','Miriam','Hera','Klymene','Artemis','Dione','Felicitas','Lydia','Ate','Iphigenia','Kassandra','Thyra','Sirona','Lomia','Althaea','Lachesis','Alkeste','Johanna','Nemesis','Antigone','Elektra','Cyrene','Sophrosyne','Hertha','Meliboea','Juewa','Siwa','Lumen','Polana','Adria','Vibilia','Adeona','Lucina','Protogeneia','Gallia','Nuwa','Bertha','Xanthippe','Aemilia','Una','Laurentia','Erigone','Eva','Loreley','Ophelia','Baucis','Ino','Phaedra','Andromache','Iduna','Irma','Klytaemnestra','Eucharis','Dejopeja','Eunike','Lamberta','Kolga','Nausikaa','Prokne','Eurykleia','Philomela','Ampella','Dynamene','Penelope','Chryseis','Pompeja','Martha','Hersilia','Hedda','Dido','Isabella','Isolda','Medea','Lilaea','Kleopatra','Eudora','Bianca','Eos','Lucia','Rosa','Oceana','Philosophia','Athamantis','Vindobona','Russia','Asterope','Carolina','Honoria','Hypatia','Vanadis','Germania','Vera','Asporina','Eukrate','Bettina','Clementina','Mathilde','Oppavia','Walpurga','Silesia','Tyche','Aletheia','Prymno','Libussa','Aline','Tirza','Adorea','Justitia','Anahita','Penthesilea','Sapientia','Adelheid','Emma','Amalia','Iclea','Nephthys','Brasilia','Felicia','Bavaria','Josephina','Olga','Nike','Polyxo','Chaldaea','Rosalia','Magdalena','Phaeo','Bamberga','Tamara','Gudrun','Svea','Etheridgea','Badenia','Roberta','Lacadiera','Devosa','Budrosa','Endymion','Desiderata','Tercidina','Hermentaria','Pariana','May','Dembowska','Ornamenta','Eleonora','Liguria','Ninina','Apollonia','Carlova','Havnia','Padua','Corduba','Vincentina','Haidea','Aeria','Bohemia','Palma','Melusina','Ursula','Campania','Huenna','Fiducia','Dodona','Ilmatar','Siegena','Aquitania','Charybdis','Industria','Wilhelmina','Lampetia','Delia','Chloe','Arsinoe','Thia','Arachne','Aspasia','Chloris','Xanthe','Elisabetha','Palatia','Vaticana','Aurelia','Diotima','Gratia','Cornelia','Hippo','Lotis','Nephele','Zeuxo','Ohio','Bathilde','Eichsfeldia','Gyptis','Valentine','Hamburga','Patientia','Mathesis','Bruchsalia','Megaira','Alekto','Lina','Argentina','Papagena','Hedwig','Tergeste','Caprera','Hansa','Genua','Venetia','Kreusa','Comacina','Veritas','Carina','Gismonda','Virtus','Tokio','Urhixidur','Evelyn','Cava','Marion','Princetonia','Iolanda','Mabella','Davida','Centesima','Armida','Amherstia','Edith','Brixia','Fidelio','Euryanthe','Turandot','Herculina','Montague','Friederike','Pamina','Deborah','Herodias','Praxedis','Ortrud','Sigelinde','Peraga','Carmen','Nanon','Suleika','Eleutheria','Cheruskia','Misa','Emanuela','Happelia','Sidonia','Klotilde','Semiramis','Bilkis','Thekla','Croatia','Irmgard','Titania','Scheila','Octavia','Luisa','Nerthus','Marianna','Tekmessa','Juvisia','Jenny','Fulvia','Valeria','Ginevra','Elfriede','Notburga','Philippina','Ute','Vundtia','Erika','Moira','Latona','Brambilla','Zelinda','Beagle','Gerlinde','Judith','Sabine','Denise','Carnegia','Rachele','Ludmilla','Melitta','Pax','Genoveva','Lanzia','Wratislavia','Lehigh','Zerbinetta','Ekard','Leonora','Galilea','Alauda','Interamnia','Erminia','Fringilla','Boliviana','Erida','Benda','Marghanna','Alagasta','Mandeville','Cantabia','Eugenisis','Aguntina','Marlu','Winchester','Faina','Sulamitis','Malabar','Lilliana','Mancunia','Massinga','Pulcova','Gedania','Tatjana','Tanete','Irmintraud','Armor','Berbericia','Theobalda','Nina','Armenia','Pickeringia','Bredichina','Hohensteina','Ani','Metcalfia','Fini','Hispania','Hormuthia','Tauris','Juliana','Adriana','Lindemannia','Burnhamia','Seraphina','Naema','Lipperta','Ara','Altona','Aida','Fatme','Lova','Rotraut','Washingtonia','Gunhild','Leopoldina','Rockefellia','Rhoda','Maritima','Palisana','Jovita','Toni','Alphonsina','Hildrun','Thuringia','Begonia','Hel','Caia','Li','Camelia','Arne','Angelica','Alsatia','Cohnia','Philippa','Aidamina','Anacostia','Gunila','Amelia','Arago','Flammario','Thomana','Vitja','Arctica','Pafuri','Amata','Asta','Feodosia','Gotho','Ljuba','Brita','Amaryllis','Nata','Freda','Lictoria','Lorraine','Sabauda','Colchis','Volga','Rusthawelia','Terentia','Aletta','Centenaria','Dysona','Pamela','Sniadeckia','Letaba','Utopia','Scythia','Spiridonia','Uzbekistania','Nyanza','Khama','Prieska','Sundmania','Mombasa','Salonta','Linzia','Imatra','Union','Strobel','Nevanlinna','Angola','Smuts','ITA','Makover','Dufour','Inge','Jiangxi','Cucula','Filipenko','Caltech','Tumilty','Kathleen'};
options.MaxConsecutiveRes = 0*ones(1, length(options.Targets)); %The maximum number of resonance orbits to each target (set to -1 to ignore)
options.MaxVisits = ones(1, length(options.Targets));           %The maximum nubmer of visists to each target (set to -1 to ignore)                    
options.AttributeIDIndex = [11 10];                             %Index of the attributes that determine the unique ID
options.RootAttrib = [0 10957.5];                                %Attributes of the root  
options.NodeCheckBoundaries = [3 2 2 2*365];               %The values used by the MyCreatedNodeCheck file. In this case, it denotes [max dV_dep, min a_per, C for the LT check, max waiting time]  
fitnessfcn = @MyCostFunctionMainBelt;                                   %The function reference to the cost function
options.NodeAttributes = @MyAttributesMainBelt;                         %The class that contains the node attributes
options.MyAttributeCalcFile = @MyAttributeCalcsMainBelt;                %The file that does the additonal calculations wrt the attributes
options.MyNodeIDCheck = @MyNodeCheckMainBelt;                           %The function that checks whether a node can be linked. Can only use the UID
options.MyCreatedNodeCheck = @MyCreatedNodeCheckMainBelt;               %After the node has been found valid using its UID and its structure has been generated, this function checks whether the node itself matches the boundaries
options.MyBestChainFile = @MyBestChainMainBelt;
options.EndTarget = {};
options.RootName = 'Start';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sets input             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofvalues = 10:10:2*365;             %Set the value of the sets. 
sets.tof = mat2cell(ones(length(options.Targets),1)... %Input should be a cell array where each line depicts a target.
   *tofvalues,[ones(length(options.Targets),1)],...    %For this mission, the ToF and the arrival epochs have been used
   [length(tofvalues)]);
load('AsteroidMainBelt/epochsnode.mat')
sets.epochsnode = epochsnode(2:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Run the optimiser        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BestSolution, BestCost, exitflag, output] = optimise_aidmap(fitnessfcn,sets,options);    
%[output] = optimise_aidmap(fitnessfcn,sets,options);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Display the result         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhysarumTreePlot(output.ListNodes)
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Display the result         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhysarumTreePlot(output.ListNodes)
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);


%Find solutions with the most asteroids
for i = 1:length(output.Solutions.Nodes);
    asteroidnum(i) = length(output.Solutions.Nodes{i});
end
maxasteroidnumindex = find(asteroidnum==max(asteroidnum));

%Plot the solutions with the most asteroids
for i = 1:length(maxasteroidnumindex)
    AllBestSolutions{i,1} = output.Solutions.Nodes{maxasteroidnumindex(i)};
    %AllBestCosts{i,1} = output.Solutions.Costs{maxasteroidnumindex(i)};
end

%Find the individual costs corresponding to the best solution
for i = 1:length(BestSolution)
    for j = 2:length(BestSolution{i})
      bestnodecosts{i}(j-1) = output.ListNodes.(char(BestSolution{i}(j))).length;
    end
end

[r] = PlotTrajectories(BestSolution,bestnodecosts,output.ListNodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Save the result          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(AllBestSolutions)
    filename = strcat(['MainBelt_',num2str(length(AllBestSolutions{1})-1),'Asteroids',num2str(i),'_',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','_',datestr(now,'yyyymmdd_HHMMSS'),'_','NewRam']);
    SaveTrajectorySolution(AllBestSolutions{i},output.ListNodes,strcat(filename));
end

%Notes:
%20160510 - GrowthFactor has been changed

