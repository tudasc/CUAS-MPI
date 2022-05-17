#include <array>

#include "PETScGrid.h"

#define NY 10
#define NX 20

std::array<std::array<PetscScalar, NX>, NY> Se = {
    {{2.94161430e+01, 2.89141085e+01, 8.97854313e+01, 1.82377004e+01, 2.43439315e+01, 4.07156863e+00, 1.33105177e+01,
      8.60147092e+01, 7.26232850e+01, 1.20097881e+01, 3.88384999e+01, 6.53053067e+00, 3.05533716e+01, 4.85085700e+01,
      2.76220506e+01, 3.49319062e+01, 4.28599208e+00, 6.58932343e+00, 8.68091891e+01, 3.88473973e+01},
     {3.87013871e+01, 1.71303342e+01, 7.02293401e+01, 5.48328630e+01, 9.21693392e+01, 6.24373725e+01, 2.34347384e+01,
      8.30452028e+00, 7.97309581e+01, 2.34022180e+01, 2.29667559e+01, 7.59988572e+01, 1.88011487e+01, 3.47992929e+01,
      4.06084125e+01, 8.68375010e+01, 1.80818441e+01, 6.58741409e+01, 9.73764997e+01, 6.58143695e+01},
     {2.88224521e+01, 2.89927575e+01, 2.92744004e+01, 2.79474205e+01, 9.85962016e+01, 3.39224482e+01, 9.29380952e+00,
      1.82402580e+01, 5.18674620e-01, 1.50059373e+01, 8.74749362e+00, 4.93290141e+01, 1.78827491e+01, 4.95980814e+01,
      8.85961401e+01, 2.36268222e+01, 6.10177578e+01, 8.01230987e+00, 4.87671510e+01, 9.53164308e+01},
     {6.34433675e+01, 3.33401424e+01, 2.46956083e+01, 2.41019150e+01, 6.02224873e+00, 3.54253086e+01, 9.85051237e+01,
      2.87027516e+01, 6.34225120e+01, 2.80644439e+01, 8.14453677e+01, 7.91444908e+01, 2.08480023e+00, 5.48941474e+01,
      5.23181066e+01, 2.72359420e+01, 1.80564799e+01, 3.50163458e+01, 8.99393212e+01, 4.84322945e+01},
     {1.85382838e+01, 7.41547383e+01, 8.81255313e+01, 9.31368054e+01, 9.72562586e+01, 5.93945607e+00, 5.98946032e+01,
      7.19974090e+01, 9.44695576e+01, 2.31192354e+01, 9.19356242e+01, 6.18920581e+01, 7.50403502e+01, 2.73075000e+01,
      4.97394697e+01, 9.60003851e+01, 2.42385667e+01, 3.90641860e+01, 5.00159851e+01, 2.32207437e+01},
     {9.06646556e+01, 8.30247230e+01, 7.70152466e+01, 3.56202445e+01, 2.84501227e+01, 7.03973525e+01, 7.91796549e+01,
      8.56232346e+01, 6.81650938e+01, 5.47531596e+01, 3.73742189e+01, 9.87938519e+01, 4.23494548e+01, 8.92535966e-02,
      7.05701077e+01, 9.21866751e+01, 1.43419544e+01, 6.99401205e+00, 3.62881036e+01, 3.30784401e+01},
     {2.66600992e+00, 7.73500613e+01, 5.23281720e+01, 3.78311200e+01, 3.91938626e+01, 2.25727183e+01, 5.89303826e+01,
      5.68292366e+01, 5.21259575e+01, 2.65173751e+01, 4.80065401e+01, 3.63201122e+01, 3.35422451e+01, 2.02990652e+01,
      7.24160254e+01, 6.51709522e+01, 2.36500856e+01, 7.77819004e+01, 5.85302067e+01, 5.86842289e+01},
     {3.73906660e+01, 8.02458098e+01, 9.74018038e+01, 7.50909165e+01, 1.07863211e+01, 5.00954441e-01, 2.63091880e+00,
      8.25646326e+01, 7.72414791e+01, 8.33135798e+00, 3.00000567e+01, 8.86667310e+01, 7.54608817e+01, 7.27668237e+01,
      3.87081547e+01, 3.19765636e+01, 5.19650990e+01, 7.83584449e+01, 4.61822429e+01, 6.56171699e+01},
     {9.21652747e+01, 5.60386170e+01, 2.69588739e+01, 9.87997135e+01, 4.60256062e+01, 8.79409216e+01, 8.40289382e+01,
      7.47223715e-01, 9.97088990e+01, 2.68689829e+01, 2.37945143e+01, 7.42336762e+01, 5.08208203e+01, 2.75234419e+01,
      7.16483520e+01, 5.44452183e+01, 9.07986483e+00, 1.62825036e+01, 9.32027239e+01, 1.61826480e+01},
     {6.66270825e+01, 8.77193382e+01, 4.26754594e+01, 1.81775998e+01, 2.17254454e+01, 8.69714876e+01, 1.68148030e+01,
      3.41058972e+01, 9.10801266e+01, 5.65776750e+01, 5.13472851e+01, 5.64316035e+01, 4.52109862e+01, 1.39374514e+01,
      7.74026393e+00, 4.56585575e+01, 2.53726849e+01, 7.28512499e+01, 7.08259854e+01, 7.21722754e+01}}};

std::array<std::array<PetscScalar, NX>, NY> TeffPow = {
    {{60.68135547, 89.98224915, 47.52978122, 73.67023354, 90.60483426, 97.49608073, 5.20037165,
      71.30247531, 19.60433394, 24.29580756, 37.60942852, 45.91205863, 38.37228151, 94.48507663,
      49.81136392, 5.59058284,  93.75887217, 86.51531015, 89.18843708, 88.60652686},
     {9.17929226,  31.39273038, 53.93505063, 93.6960327,  52.82987497, 32.0527514,  73.42238915,
      46.71278264, 20.65822352, 1.68386586,  4.95952478,  22.37292686, 48.45713146, 62.145702,
      22.29040969, 89.25479229, 52.61110359, 22.65352636, 48.59269249, 66.03999027},
     {61.23842332, 71.88582593, 19.58645317, 97.33072769, 11.036336,   87.09710827, 36.25284102,
      68.46491615, 63.0743987,  39.45661385, 36.1720578,  83.73411868, 46.89061859, 84.69869315,
      70.52764993, 17.44368911, 83.08772884, 12.77070925, 13.36103111, 78.54710271},
     {96.44942334, 61.70070712, 25.92356087, 6.33366667,  72.84720738, 77.43038589, 80.30811883,
      7.67297831,  6.32457338,  19.46688868, 42.57031725, 25.12676205, 88.92646342, 86.73140024,
      21.81475927, 37.55953289, 28.570403,   9.56321525,  10.21877706, 8.5801699},
     {74.74565108, 65.97710098, 34.36648359, 52.02449694, 26.24473474, 94.8344978,  67.17641608,
      12.94304264, 58.93947433, 3.10025945,  89.14260007, 25.94078815, 60.40622692, 1.93870785,
      76.2056459,  1.97754578,  83.11012923, 20.62213525, 79.21517918, 68.88763808},
     {47.84308753, 37.45553281, 30.48834879, 76.38495592, 35.98447911, 90.61480182, 91.70508372,
      89.19971556, 51.49557168, 18.38458375, 9.99231195,  69.6337708,  49.75562359, 14.276455,
      90.66851169, 55.05299308, 91.60370858, 4.19417423,  11.26718896, 33.60320595},
     {25.82726535, 88.77852959, 53.8286781,  15.38450237, 71.08492098, 56.57979208, 96.09665016,
      79.21865204, 88.7647144,  80.07106664, 30.83456706, 99.69444968, 53.93856345, 25.57358272,
      24.84531085, 70.43149798, 60.60739068, 85.66940796, 57.18117295, 67.68890213},
     {17.57516755, 8.23445567,  89.20626904, 39.863904,   9.69156607,  15.96900742, 31.42712893,
      46.58305481, 40.71355373, 20.03367259, 95.13993145, 53.67586958, 58.76657393, 35.08710936,
      12.15481832, 69.68244688, 26.06816414, 10.43666554, 27.79405803, 92.1911942},
     {65.25065033, 26.84303477, 22.79122176, 9.7903579,   47.95502493, 66.58074533, 86.64225326,
      75.17105102, 47.62372316, 38.90584528, 75.12833522, 94.81118643, 51.29925369, 88.79491456,
      85.29938626, 98.89257693, 78.19386205, 94.61121583, 41.66165042, 76.51775988},
     {0.41777659,  19.06178583, 83.11287479, 15.76921537, 38.02785595, 11.2274433,  61.94517324,
      84.0354684,  45.50979813, 36.20947539, 58.49475378, 82.67780511, 73.37257522, 30.37875836,
      77.74269443, 51.39481683, 17.73788486, 25.03282614, 3.40935245,  58.67038533}}};

std::array<std::array<PetscScalar, NX>, NY> u_n = {
    {{28.09387556, 66.34389685, 43.98926676, 55.15831325, 57.83683051, 10.45031262, 41.24323809,
      68.80014622, 12.85353449, 2.0156961,   16.97804254, 95.5785155,  53.83963516, 21.53842467,
      14.81691747, 51.53104372, 55.64217153, 20.23872146, 86.01555489, 19.32046073},
     {28.92829482, 5.81556337,  5.27494935, 16.07288859, 35.29536806, 49.02446322, 17.88484884,
      11.17914247, 42.66914442, 7.37966386, 35.85269581, 69.57705409, 49.19151783, 53.70745627,
      15.28393447, 76.25023623, 2.10792218, 1.63631078,  77.48081464, 84.4795664},
     {0.62860842,  9.05840173,  12.7286033,  0.52601021,  39.39981882, 84.64877837, 16.7947071,
      16.30071632, 27.37622894, 90.36696577, 67.15237107, 96.76657744, 38.43661436, 56.28684077,
      78.04449169, 70.26825629, 22.18454319, 7.4480884,   93.27187813, 95.9204784},
     {12.4125356,  75.0543057,  1.87608478,  67.90266451, 63.66869108, 30.72702625, 20.44184374,
      83.97502821, 98.89276471, 8.27723966,  34.44743608, 20.49032975, 69.95914646, 45.63752997,
      93.70826077, 77.99684146, 86.96436919, 9.24589587,  24.7966672,  64.75132488},
     {62.40419063, 74.59303786, 2.59262374,  16.15427984, 77.7258353,  51.86866523, 15.7866676,
      47.44818248, 13.66772346, 78.5304741,  83.81728781, 52.56117198, 48.36148179, 70.16861946,
      97.94012039, 98.20682273, 10.38228253, 72.88895276, 92.18249289, 16.72255978},
     {79.25184062, 36.03739872, 54.91632443, 15.31594857, 87.23753033, 33.15896144, 77.38891488,
      7.8537693,   57.76081853, 16.27221065, 71.29728086, 77.86673342, 74.68815583, 9.07744639,
      69.54490626, 18.94560504, 36.96767777, 90.32428701, 89.37997193, 35.54253641},
     {51.18224948, 31.70104721, 28.28405059, 31.88240105, 31.81022349, 87.0505811,  15.51908143,
      23.22102871, 74.73881663, 3.03687392,  3.55463826,  45.36457623, 97.38424118, 12.56854422,
      55.29083783, 16.25937373, 82.05734698, 3.12202079,  69.61723527, 49.39940728},
     {73.9302183,  82.35852055, 17.37425406, 99.6212366,  69.51704233, 64.1092411,  96.56071817,
      19.47350029, 6.92748155,  99.82432492, 59.81191679, 25.45557997, 86.95813419, 98.67229634,
      68.52679166, 97.26718314, 42.11219036, 38.1830747,  77.81208586, 62.45191734},
     {76.10424052, 56.28355187, 98.30589259, 61.09311947, 32.37122045, 62.60610944, 22.05682466,
      18.09815215, 21.46346098, 84.21059768, 0.24036326,  80.10659927, 45.43920055, 66.34209541,
      76.75051261, 20.61492172, 26.523529,   38.62224297, 71.48987175, 91.80872048},
     {56.08902535, 87.96283575, 63.19303124, 35.41950103, 50.37683678, 69.81383351, 52.73571757,
      50.81591431, 91.66922317, 51.70122995, 26.8727926,  56.04088529, 93.7225918,  15.36705651,
      74.27585625, 82.28208885, 54.96862011, 28.94664115, 57.06853161, 67.5220487}}};

std::array<std::array<PetscScalar, NX>, NY> current_Q = {
    {{11.80319311, 41.07773885, 93.05962722, 78.92108556, 48.9113483,  97.49591666, 30.6395652,
      15.71459821, 14.56212443, 46.78484117, 79.46024987, 30.45826421, 74.1408329,  7.80921743,
      80.5448044,  75.61297889, 81.98163269, 16.85658392, 93.21797538, 43.02412003},
     {99.20625059, 9.29062304,  17.06772982, 63.95993231, 12.57783805, 84.02830845, 90.6657825,
      77.96345521, 95.14118294, 8.70976247,  35.64036563, 94.80801997, 97.27281221, 75.6631087,
      11.55866615, 1.38543396,  89.88497011, 50.05200106, 44.19747494, 82.08217832},
     {73.90456184, 24.00271479, 11.87145952, 92.48117093, 56.72476213, 44.18638903, 83.27694886,
      16.00293687, 28.68114068, 97.39469783, 35.17730108, 94.25010986, 38.09907479, 57.28556251,
      75.70129835, 12.16437951, 58.46195315, 67.07392158, 20.69033463, 47.50282809},
     {86.72200915, 90.02476045, 82.86212348, 35.20416334, 93.6256578,  67.64342693, 37.78010029,
      12.44886509, 35.55769941, 41.95021341, 26.61329071, 91.66574571, 82.50725209, 69.98483551,
      8.68195616,  0.44238022,  88.23291845, 60.43460639, 54.08850493, 34.18910395},
     {17.7256557,  42.59646212, 99.96678577, 36.51826181, 6.68411993,  21.2993519,  52.6342458,
      35.31474229, 52.11533305, 35.66713065, 34.19832161, 89.31893299, 23.45860365, 51.88735734,
      74.48325824, 95.14746264, 83.19686001, 53.71675708, 93.08345332, 70.57607281},
     {55.52741831, 73.72903607, 8.98026261,  99.25460244, 24.38564587, 3.47024342,  64.62957869,
      14.9629595,  70.18179547, 17.90880797, 48.35716608, 42.42309694, 77.91316768, 85.4766328,
      86.41111209, 46.87126851, 4.0131813,   99.51129523, 62.0332488,  55.99651438},
     {32.35639075, 65.09257045, 40.55067875, 76.15776058, 75.10605439, 12.61893246, 52.95483312,
      31.77578541, 19.25432321, 26.18721578, 68.71682008, 88.56128645, 14.41891219, 59.42567364,
      95.48236088, 84.00636406, 58.24372768, 75.08825975, 32.32471593, 48.64289568},
     {69.66349679, 50.16375037, 68.62882208, 78.06146608, 36.52544343, 12.62404649, 29.25223319,
      40.72940505, 51.99382108, 74.14090265, 67.34417334, 56.18106125, 8.06213214,  56.65093786,
      20.33595884, 20.07858153, 26.40641418, 33.73820201, 56.27879561, 94.01190966},
     {61.99547479, 71.27756592, 93.65210188, 40.61315691, 17.09516822, 69.73801531, 45.27220814,
      5.8764625,   25.0923235,  91.82511527, 56.97886316, 71.09192002, 97.51582851, 51.25113983,
      38.52819752, 79.2530608,  22.79008734, 92.30300215, 60.19928455, 3.15245726},
     {68.29616684, 6.51455461,  9.5650454,   67.65454825, 82.17471523, 47.87206249, 36.25795578,
      12.24954396, 23.65935172, 58.24966295, 35.25689647, 82.5388064,  13.39163619, 0.17364557,
      46.18177003, 60.60350219, 6.27916636,  18.35530765, 53.26105516, 58.66903733}}};

std::array<std::array<PetscScalar, NX>, NY> dirich_val = {
    {{8.29842071e+00, 5.66585089e+00, 1.91771991e+00, 2.51873477e+00, 6.75669113e+00, 8.02822026e+00, 7.11489931e+00,
      6.21650194e-01, 7.92404011e+00, 7.03534198e+00, 1.85145803e+00, 2.78505066e+00, 5.91234450e-01, 3.64215513e+00,
      5.53774278e+00, 5.06391010e+00, 1.94474974e+00, 5.43735539e+00, 8.89543147e+00, 4.94238867e+00},
     {1.41810062e-01, 8.19961889e+00, 9.06026724e+00, 9.25102103e+00, 6.66746716e+00, 6.13044126e+00, 6.66198869e+00,
      1.40502003e+00, 1.21035803e+00, 5.59510319e+00, 6.30410637e+00, 6.35655481e+00, 2.39462229e+00, 1.03974595e+00,
      8.06368950e+00, 7.04016690e+00, 3.86125100e-01, 9.51978056e-01, 5.44699054e+00, 9.02462881e+00},
     {4.43698010e+00, 8.32003102e+00, 8.29631276e+00, 1.69550488e+00, 3.79841664e+00, 3.12158160e+00, 8.54379293e+00,
      2.85892579e+00, 1.08679361e-02, 4.99047777e+00, 1.28075440e+00, 3.63250312e+00, 8.58864857e-02, 6.98424018e+00,
      7.67802510e+00, 4.32226396e+00, 4.62065153e+00, 5.37203855e+00, 2.65782963e+00, 6.77843541e-01},
     {7.28768120e+00, 5.97401644e+00, 5.83804845e-02, 6.30291932e+00, 3.30424883e+00, 8.39785573e+00, 1.38179220e+00,
      4.36020442e+00, 4.82776991e+00, 8.92879864e+00, 7.54992161e+00, 5.31214730e-01, 6.58666525e+00, 6.66053824e+00,
      2.75410043e+00, 4.48381684e+00, 2.60294575e+00, 7.81603555e+00, 5.16386742e+00, 8.13159913e+00},
     {6.51129860e+00, 7.88284573e+00, 1.80746481e+00, 6.74327841e+00, 7.96838274e-01, 7.26717986e+00, 6.16579848e+00,
      8.46476331e+00, 6.56876419e+00, 3.48119591e+00, 4.12908147e+00, 3.43621714e+00, 9.61954610e+00, 6.01949659e+00,
      4.33947884e+00, 2.76523389e+00, 8.69674829e+00, 5.57836293e+00, 5.40886878e+00, 4.40192989e+00},
     {6.94494799e-01, 4.64273695e+00, 2.17906679e+00, 5.29358051e-01, 9.31229076e+00, 4.13983341e+00, 2.01378870e+00,
      1.00441023e-01, 7.24526421e+00, 8.21420150e+00, 3.32772086e+00, 5.87673277e+00, 9.50267622e+00, 9.79303376e+00,
      8.67294360e+00, 4.83769361e+00, 9.17130258e+00, 5.55335379e-01, 3.48044261e+00, 2.62638616e+00},
     {5.85061679e+00, 3.92574210e+00, 9.28571166e+00, 4.12552502e+00, 4.07411024e+00, 8.11551752e+00, 6.55869671e+00,
      8.98803980e+00, 1.57877377e+00, 4.08172566e-01, 2.72285240e+00, 1.15847259e+00, 1.51950674e+00, 8.88664465e+00,
      5.40741747e+00, 2.06799632e+00, 5.98809196e+00, 8.04467534e+00, 6.92511416e+00, 8.74507150e+00},
     {9.40183491e+00, 9.01449528e+00, 4.82927201e+00, 2.08333428e+00, 6.86158714e+00, 4.88020328e+00, 2.39376150e+00,
      7.41591009e+00, 7.31334714e+00, 9.84054676e+00, 4.19099827e+00, 7.45292514e+00, 8.52844375e+00, 3.84588589e+00,
      7.03331596e+00, 5.73341265e+00, 3.66553696e+00, 9.59226880e+00, 2.90485550e+00, 5.52029348e+00},
     {4.43673318e+00, 7.78943260e+00, 2.43361211e+00, 5.90077746e+00, 1.37997403e+00, 6.92213714e+00, 7.11428907e+00,
      1.65254575e-01, 7.69932750e+00, 8.68042420e+00, 9.72002740e+00, 7.24116339e-01, 8.70489820e+00, 6.55875259e+00,
      5.55889467e+00, 5.09264649e+00, 4.81189145e+00, 4.94203908e+00, 8.14382484e+00, 1.05316410e+00},
     {4.83297115e+00, 5.25562669e-01, 8.94005528e+00, 7.97932308e-01, 7.17189076e+00, 3.28779717e+00, 7.06498343e-03,
      8.76249847e+00, 5.65596981e+00, 3.50448752e+00, 6.25347474e+00, 7.94681317e+00, 6.98720844e+00, 5.29540758e+00,
      7.04815711e+00, 8.06750562e+00, 6.60220556e+00, 7.84214889e-01, 3.07540100e+00, 2.00089401e-01}}};

std::array<std::array<PetscScalar, NX>, NY> dirich_mask = {
    {{true, true, true, true, true, true, true, true, true, true,
      true, true, true, true, true, true, true, true, true, true},
     {true,  false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, true},
     {true,  false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, true},
     {true,  false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, true},
     {true,  false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, true},
     {true,  false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, true},
     {true,  false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, true},
     {true,  false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, true},
     {true,  false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, true},
     {true, true, true, true, true, true, true, true, true, true,
      true, true, true, true, true, true, true, true, true, true}}};

std::array<std::array<PetscScalar, NX>, NY> bnd_mask = {
    {{2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0},
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0},
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0},
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0},
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0},
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0},
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0},
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0},
     {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0},
     {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}}};

std::array<PetscScalar, NX *NY> b = {
    8.29842071e+00,  1.41810062e-01,  4.43698010e+00,  7.28768120e+00,  6.51129860e+00,  6.94494799e-01,
    5.85061679e+00,  9.40183491e+00,  4.43673318e+00,  4.83297115e+00,  5.66585089e+00,  -1.31732249e+05,
    -7.21226264e+04, 3.93888845e+05,  8.75337944e+04,  -3.62937915e+03, 1.78634191e+04,  4.70087007e+04,
    1.06819847e+04,  5.25562669e-01,  1.91771991e+00,  -2.38273688e+04, 4.29635071e+04,  1.74861416e+04,
    -4.44329215e+03, 5.79691258e+04,  2.17064166e+04,  -4.75485911e+04, 3.11538641e+05,  8.94005528e+00,
    2.51873477e+00,  1.32222682e+04,  3.37169344e+04,  1.29328843e+05,  -2.52163954e+03, -1.71233091e+04,
    5.99994638e+04,  1.18256648e+05,  1.38607525e+04,  7.97932308e-01,  6.75669113e+00,  -1.70882548e+03,
    1.90778154e+04,  9.40104193e+05,  3.15455244e+04,  2.96278574e+05,  -6.20121101e+04, 1.99044882e+05,
    -4.65365137e+04, 7.17189076e+00,  8.02822026e+00,  8.17231853e+04,  3.00710787e+05,  -1.02359095e+05,
    8.64809861e+05,  -1.38796630e+05, 4.83006801e+05,  -1.05154483e+06, 7.37903781e+04,  3.28779717e+00,
    7.11489931e+00,  1.15262980e+05,  4.14739627e+04,  5.83404607e+03,  -8.08270122e+04, 2.45625659e+05,
    -1.74421325e+05, 3.77175607e+06,  -5.86840339e+04, 7.06498343e-03,  6.21650194e-01,  -8.95914325e+04,
    -2.64785286e+04, 7.57133689e+04,  4.15920481e+04,  -1.03327647e+05, -7.71705442e+03, -4.39487992e+03,
    -3.12266799e+06, 8.76249847e+00,  7.92404011e+00,  6.91992768e+04,  -3.40871937e+06, 5.77915111e+04,
    -1.61071243e+04, 1.08887770e+05,  2.50504112e+05,  -5.07160562e+04, -3.12108826e+04, 5.65596981e+00,
    7.03534198e+00,  1.06115737e+02,  6.01983149e+05,  -5.69038505e+04, 9.90295955e+04,  -1.28676109e+04,
    -2.50567627e+05, 1.06764263e+06,  3.71535677e+05,  3.50448752e+00,  1.85145803e+00,  5.16721417e+04,
    8.65144967e+04,  -1.63155590e+04, 5.32485670e+04,  8.22554656e+04,  -3.61246719e+04, 3.22886751e+05,
    -5.04579243e+05, 6.25347474e+00,  2.78505066e+00,  4.13647613e+04,  2.22576492e+05,  -1.58898733e+04,
    4.20201791e+04,  4.93443480e+04,  3.28725160e+04,  -7.38821759e+04, 2.00072316e+05,  7.94681317e+00,
    5.91234450e-01,  1.84596830e+05,  -2.30930336e+05, 4.51660826e+06,  -2.19773875e+04, 1.03336040e+05,
    2.17842492e+05,  5.37510046e+04,  -6.30631380e+04, 6.98720844e+00,  3.64215513e+00,  1.74058926e+05,
    5.47059492e+04,  -1.92969801e+04, 8.48188244e+04,  9.51258420e+06,  -1.40140269e+05, 8.93820526e+04,
    1.08289862e+05,  5.29540758e+00,  5.53774278e+00,  -9.85604584e+04, 6.62823671e+04,  4.50838835e+04,
    1.02504352e+05,  8.54906836e+04,  7.25020889e+04,  -1.62343477e+03, 9.39585025e+04,  7.04815711e+00,
    5.06391010e+00,  7.04580702e+04,  4.84330464e+04,  -1.75013288e+04, 4.60504537e+04,  -2.18178865e+04,
    -9.09739973e+04, 3.93936303e+05,  -1.79496056e+05, 8.06750562e+00,  1.94474974e+00,  -2.36153536e+05,
    5.32866633e+03,  5.27742488e+05,  -1.26892200e+05, 2.27162304e+03,  5.97803516e+05,  -2.41490536e+04,
    -7.20541538e+04, 6.60220556e+00,  5.43735539e+00,  -7.56644036e+03, 1.99270797e+05,  2.16902690e+04,
    9.96757172e+04,  7.17444510e+05,  -8.03601135e+04, 1.72176419e+04,  2.04911125e+05,  7.84214889e-01,
    8.89543147e+00,  2.73260419e+04,  5.32624563e+04,  4.25727312e+03,  1.95320491e+05,  9.60155009e+04,
    9.29925738e+04,  8.01189808e+04,  3.20098188e+04,  3.07540100e+00,  4.94238867e+00,  9.02462881e+00,
    6.77843541e-01,  8.13159913e+00,  4.40192989e+00,  2.62638616e+00,  8.74507150e+00,  5.52029348e+00,
    1.05316410e+00,  2.00089401e-01};

std::array<PetscScalar, 776> ANonZeroPython = {
    1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,
    1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  -1.00199358e+03,
    -3.28329629e+03, 1.01682657e+04,  -3.08258819e+03, -2.79938764e+03, -2.75638400e+03, -1.82134334e+03,
    8.62934546e+03,  -2.76758109e+03, -1.28303703e+03, -2.72754497e+03, -2.40670260e+03, 8.76951159e+03,
    -2.31110485e+03, -1.32315917e+03, -1.14207619e+03, -1.03907810e+03, 3.69719439e+03,  -7.78631313e+02,
    -7.36408776e+02, -6.11511859e+02, -6.95445877e+02, 2.56394501e+03,  -7.66759350e+02, -4.89227928e+02,
    -6.25083852e+02, -8.23011405e+02, 2.73151118e+03,  -2.35435235e+02, -1.04698069e+03, -1.68868827e+02,
    -2.26939325e+02, 8.13613032e+02,  -1.89773155e+02, -2.27031725e+02, -8.20197319e+02, -2.71750112e+02,
    2.10519787e+03,  -4.80693936e+02, -5.31556506e+02, 1.00000000e+00,  1.00000000e+00,  -6.82826376e+02,
    -8.69404584e+02, 3.22558832e+03,  -4.94439396e+02, -1.17791797e+03, -1.27069320e+03, -1.18616102e+03,
    4.72492789e+03,  -9.21031796e+02, -1.34604188e+03, -1.78632227e+03, -1.09179953e+03, 4.82327919e+03,
    -1.44605086e+03, -4.98106535e+02, -6.19663783e+02, -4.05229961e+02, 2.03646901e+03,  -4.43041143e+02,
    -5.67534128e+02, -5.27402236e+02, -5.06954634e+02, 2.32990046e+03,  -6.10765886e+02, -6.83777703e+02,
    -1.54761799e+03, -8.98909393e+02, 4.55052741e+03,  -1.55042412e+03, -5.52575907e+02, -1.87043196e+02,
    -8.32950280e+02, 2.15500131e+03,  -4.50408770e+02, -6.83599064e+02, -1.10493085e+03, -1.62731673e+03,
    4.95056219e+03,  -1.60339510e+03, -6.13919512e+02, 1.00000000e+00,  1.00000000e+00,  -1.50866464e+03,
    -1.81770561e+03, 6.92031282e+03,  -2.10404579e+03, -1.48889678e+03, -1.40995370e+03, -4.12813964e+03,
    6.91046448e+03,  -5.14224751e+02, -8.57146382e+02, -5.10376203e+02, -5.96270269e+02, 2.25806575e+03,
    -5.66146819e+02, -5.84272459e+02, -5.36997660e+02, -1.46507307e+02, 1.94015085e+03,  -8.03001833e+02,
    -4.52644051e+02, -1.47840951e+03, -2.09962133e+03, 6.10739781e+03,  -8.68791110e+02, -1.65957585e+03,
    -7.64325430e+02, -8.18018387e+02, 3.10037353e+03,  -7.09109778e+02, -8.07919934e+02, -8.86708873e+02,
    -3.57252493e+02, 1.74883005e+03,  -2.52961022e+02, -2.50907666e+02, -1.67516465e+02, -1.92258402e+02,
    7.07396216e+02,  -1.47746947e+02, -1.98874402e+02, 1.00000000e+00,  1.00000000e+00,  -8.85766067e+02,
    -8.75002653e+02, 2.52420542e+03,  -2.39367960e+02, -5.23068736e+02, -2.42960986e+02, -2.23765078e+02,
    9.42734143e+02,  -2.34920014e+02, -2.40088065e+02, -2.33834333e+03, -3.84610834e+03, 2.89901549e+04,
    -7.74242869e+03, -1.50622745e+04, -4.33471547e+02, -4.79422425e+02, 1.80178879e+03,  -3.77107077e+02,
    -5.10787742e+02, -2.07782926e+03, -1.28913410e+03, 7.58519147e+03,  -2.02937515e+03, -2.18785295e+03,
    -7.79829136e+02, -1.47308707e+03, 4.72233315e+03,  -5.25881104e+02, -1.94253584e+03, -1.74673889e+03,
    -1.91087503e+03, 6.81625565e+03,  -1.80634314e+03, -1.35129859e+03, -4.26908748e+02, -4.23325161e+02,
    3.42860175e+03,  -1.11363687e+03, -1.46373098e+03, 1.00000000e+00,  1.00000000e+00,  -7.72148119e+02,
    -9.33670204e+02, 3.47730994e+03,  -9.06880461e+02, -8.63611155e+02, -6.97820250e+02, -1.66919654e+03,
    7.11182872e+03,  -2.92017337e+03, -1.82363856e+03, -2.56056382e+03, -2.79628982e+03, 1.09550972e+04,
    -2.90794724e+03, -2.68929635e+03, -8.36394849e+03, -1.73441687e+04, 6.05632026e+04,  -1.88544023e+04,
    -1.59996831e+04, -8.84190709e+02, -1.59075434e+03, 5.23633534e+03,  -1.19572540e+03, -1.56466489e+03,
    -3.37289829e+03, -3.72909905e+03, 1.22490497e+04,  -1.33335049e+03, -3.81270191e+03, -2.90955412e+04,
    -6.00800044e+04, 2.02392032e+05,  -6.21343436e+04, -5.10811429e+04, -7.66072315e+02, -3.53947568e+02,
    2.41966121e+03,  -2.64017280e+02, -1.03462405e+03, 1.00000000e+00,  1.00000000e+00,  -2.30092653e+03,
    -5.00810750e+02, 8.24961165e+03,  -2.50277059e+03, -2.94410378e+03, -6.65628928e+03, -6.31084312e+03,
    2.56263783e+04,  -6.49491938e+03, -6.16332651e+03, -9.67149216e+02, -6.12785826e+02, 2.65017161e+03,
    -8.97408393e+02, -1.71828175e+02, -1.58661064e+03, -1.47591469e+03, 5.06587087e+03,  -1.56447406e+03,
    -4.37871482e+02, -1.39111829e+03, -1.18342967e+03, 5.38786660e+03,  -1.43221298e+03, -1.38010566e+03,
    -1.46041892e+03, -1.92434063e+03, 6.13766439e+03,  -9.71184335e+02, -1.78072050e+03, -9.72638356e+03,
    -2.17537175e+04, 6.99033905e+04,  -2.11840434e+04, -1.72382460e+04, -1.08279117e+03, -6.63265528e+02,
    3.94348331e+03,  -1.03883224e+03, -1.15759437e+03, 1.00000000e+00,  1.00000000e+00,  -8.30804184e+03,
    -8.21309381e+03, 2.87709814e+04,  -8.08053321e+03, -4.16831254e+03, -3.14034936e+03, -3.67894752e+03,
    1.20840957e+04,  -9.14153990e+02, -4.34964482e+03, -5.89698013e+02, -5.80934011e+02, 1.86912877e+03,
    -4.05593623e+02, -2.91903121e+02, -3.64265034e+02, -1.61695444e+02, 1.26258080e+03,  -3.79397575e+02,
    -3.56222744e+02, -1.27624575e+03, -3.19021379e+02, 3.70194809e+03,  -1.18421070e+03, -9.21470254e+02,
    -1.84655904e+03, -1.78422159e+03, 6.65932894e+03,  -1.24743006e+03, -1.78011825e+03, -5.49296038e+02,
    -8.58606111e+02, 2.88663236e+03,  -8.41821604e+02, -6.35908603e+02, -1.30177113e+05, -9.30172451e+04,
    4.45812540e+05,  -1.28327845e+05, -9.42893373e+04, 1.00000000e+00,  1.00000000e+00,  -4.34158034e+02,
    -3.04884868e+02, 1.25891072e+03,  -4.71675766e+02, -4.71920499e+01, -1.52964191e+05, -7.25062675e+04,
    3.65348818e+05,  -2.67827685e+04, -1.13094591e+05, -1.32104871e+02, -2.19031726e+02, 7.51674671e+02,
    -2.17640240e+02, -1.81897833e+02, -2.71485494e+02, -1.46113638e+02, 1.19701202e+03,  -7.03066640e+02,
    -7.53462493e+01, -1.15747312e+03, -9.74375458e+02, 3.76856517e+03,  -1.15540072e+03, -4.80315867e+02,
    -1.94073675e+03, -1.51091706e+03, 6.69840932e+03,  -1.29403968e+03, -1.95171583e+03, -6.79732712e+02,
    -8.73275058e+02, 2.66082856e+03,  -6.86731415e+02, -4.20089374e+02, -7.06609236e+02, -5.31990132e+02,
    2.32262607e+03,  -5.64036822e+02, -5.18989880e+02, 1.00000000e+00,  1.00000000e+00,  -1.60782510e+02,
    -1.62617496e+02, 6.20984339e+02,  -1.66770825e+02, -1.29813508e+02, -3.90907231e+03, -2.60084201e+02,
    9.30872461e+03,  -2.09934504e+03, -3.03922306e+03, -4.11068808e+02, -1.12251075e+03, 2.91517988e+03,
    -2.30292933e+02, -1.15030739e+03, -3.07878990e+02, -2.79552632e+02, 1.17892311e+03,  -2.77309873e+02,
    -3.13181616e+02, -5.97970535e+02, -1.17092644e+02, 1.66173042e+03,  -6.59931307e+02, -2.85735935e+02,
    -3.83654325e+03, -1.36262824e+03, 8.68940907e+03,  -1.46039441e+03, -2.02884317e+03, -3.89472216e+03,
    -4.64820099e+03, 1.71802265e+04,  -3.83594348e+03, -4.80035992e+03, -1.92593481e+03, -1.18942419e+03,
    7.10863671e+03,  -1.68685007e+03, -2.30542764e+03, 1.00000000e+00,  1.00000000e+00,  -1.32274842e+02,
    -4.61066235e+02, 1.48045558e+03,  -4.58941524e+02, -4.27172976e+02, -5.21365234e+03, -1.20496206e+03,
    1.88008890e+04,  -5.40265338e+03, -6.97862124e+03, -3.96372907e+02, -5.80262296e+02, 2.30137939e+03,
    -8.54902474e+02, -4.68841709e+02, -7.87564077e+01, -7.57354366e+02, 1.60149112e+03,  -2.36189167e+02,
    -5.28191178e+02, -4.18602601e+02, -5.80994042e+02, 2.05361752e+03,  -4.87982985e+02, -5.65037887e+02,
    -1.12067221e+03, -3.79906214e+02, 3.85942483e+03,  -1.17229467e+03, -1.18555173e+03, -1.33311471e+03,
    -1.87592350e+03, 9.35602874e+03,  -3.38166013e+03, -2.76433040e+03, -2.60330995e+03, -4.26358758e+03,
    1.44652414e+04,  -3.34026928e+03, -4.25707458e+03, 1.00000000e+00,  1.00000000e+00,  -1.29091118e+02,
    -4.78339707e+02, 1.65657289e+03,  -5.61426774e+02, -4.86715291e+02, -1.23751601e+03, -8.64963430e+02,
    4.52291796e+03,  -9.46855811e+02, -1.47258271e+03, -4.82471806e+02, -5.90154326e+02, 2.06158103e+03,
    -3.89738336e+02, -5.98216564e+02, -7.84585084e+02, -4.98378033e+02, 2.73054483e+03,  -7.37980694e+02,
    -7.08601015e+02, -2.13756719e+02, -4.62327798e+02, 2.38984972e+03,  -1.00288418e+03, -7.09881024e+02,
    -1.56701710e+03, -2.72793185e+03, 8.94645041e+03,  -2.32156429e+03, -2.32893718e+03, -9.35300846e+02,
    -9.50970836e+02, 3.58601016e+03,  -9.34132590e+02, -7.64605886e+02, -1.36454271e+03, -1.11575349e+03,
    5.00278182e+03,  -1.43778940e+03, -1.08369622e+03, 1.00000000e+00,  1.00000000e+00,  -1.96742265e+03,
    -2.75260073e+03, 1.12839237e+04,  -3.06314797e+03, -3.49975231e+03, -4.06207418e+03, -3.22046125e+03,
    1.55113222e+04,  -4.14903837e+03, -4.07874843e+03, -2.27098716e+04, -3.55891233e+04, 1.50895007e+05,
    -4.16977425e+04, -5.08972695e+04, -5.84442572e+02, -1.15846292e+03, 2.68305228e+03,  -8.78652110e+02,
    -6.04946822e+01, -1.65602795e+03, -1.55691171e+03, 5.32391220e+03,  -1.47692507e+03, -6.33047471e+02,
    -2.52181270e+03, -1.86472227e+03, 7.66380496e+03,  -2.02634344e+03, -1.24992655e+03, -8.98413892e+02,
    -9.00706522e+02, 3.38089219e+03,  -8.77174884e+02, -7.03596894e+02, -1.58294876e+03, -1.30246599e+03,
    5.86825010e+03,  -1.43566778e+03, -1.54616757e+03, 1.00000000e+00,  1.00000000e+00,  -1.89082473e+03,
    -2.60343021e+03, 8.12389627e+03,  -2.48931102e+03, -1.13933031e+03, -1.47060598e+03, -1.74656480e+03,
    7.18122395e+03,  -2.08795087e+03, -1.87510230e+03, -1.93300459e+03, -1.88650999e+03, 4.67136991e+03,
    -8.34842582e+01, -7.67371072e+02, -1.66237925e+02, -1.67821924e+02, 6.53436190e+02,  -1.51059707e+02,
    -1.67316634e+02, -3.00371260e+05, -4.62173302e+04, 9.28631894e+05,  -2.48071367e+05, -3.33970937e+05,
    -2.06538293e+03, -1.09075278e+03, 6.41852538e+03,  -1.76106308e+03, -1.50032659e+03, -7.29646277e+02,
    -4.91266933e+02, 2.35697231e+03,  -8.35243533e+02, -2.99815568e+02, -2.85493015e+03, -2.20822741e+03,
    1.08716202e+04,  -1.98744076e+03, -3.82002188e+03, 1.00000000e+00,  1.00000000e+00,  -9.76346690e+02,
    -9.16439229e+02, 3.96321270e+03,  -1.00797118e+03, -1.06145561e+03, -1.04972380e+03, -4.62007818e+02,
    2.34868018e+03,  -4.54478123e+02, -3.81470444e+02, -8.05154918e+02, -7.69618973e+02, 2.99662223e+03,
    -7.83408373e+02, -6.37439969e+02, -9.18586187e+01, -8.24022512e+02, 3.02227662e+03,  -2.01174322e+03,
    -9.36522670e+01, -4.22389993e+02, -1.41792388e+03, 3.68218113e+03,  -6.67829091e+02, -1.17303817e+03,
    -4.20559222e+02, -6.50805821e+02, 1.95767117e+03,  -2.72379318e+02, -6.12926813e+02, -5.63618358e+02,
    -5.09572925e+02, 2.38456511e+03,  -6.64217749e+02, -6.46156080e+02, -1.46744688e+03, -3.58844867e+02,
    4.74390742e+03,  -1.37188822e+03, -1.54472746e+03, 1.00000000e+00,  1.00000000e+00,  -4.96375721e+02,
    -1.46414572e+02, 1.97105775e+03,  -4.06091102e+02, -9.21176353e+02, -1.43044243e+03, -1.49253828e+03,
    5.61700903e+03,  -1.21838411e+03, -1.47464422e+03, -1.22447214e+03, -1.05693222e+03, 3.88894113e+03,
    -1.66693589e+02, -1.43984318e+03, -4.85228688e+01, -4.72920698e+01, 1.93495824e+02,  -4.80558405e+01,
    -4.86250452e+01, -8.97976086e+02, -5.00438831e+01, 2.66051730e+03,  -8.10044361e+02, -9.01452971e+02,
    -6.81066060e+02, -1.14583712e+03, 4.33477162e+03,  -1.29889498e+03, -1.20797346e+03, -7.82182533e+02,
    -2.64725828e+03, 7.95366162e+03,  -3.08945031e+03, -1.43377051e+03, -2.03281720e+03, -1.81448449e+03,
    7.28769007e+03,  -1.50113540e+03, -1.93825297e+03, 1.00000000e+00,  1.00000000e+00,  -4.42392116e+03,
    -4.50416419e+03, 1.53509011e+04,  -4.30541100e+03, -2.11640474e+03, -5.71000279e+02, -1.27585433e+03,
    3.12829523e+03,  -8.42026577e+02, -4.38414046e+02, -2.17182338e+03, -2.84543687e+03, 8.82284862e+03,
    -2.84563312e+03, -9.58955247e+02, -1.92586596e+02, -2.11984966e+03, 8.30540774e+03,  -4.34462126e+03,
    -1.64735023e+03, -5.79432551e+03, -7.34261101e+03, 1.99598803e+04,  -6.14614967e+03, -6.75794143e+02,
    -3.32873132e+03, -3.72716614e+03, 1.25466476e+04,  -1.86262674e+03, -3.62712344e+03, -8.82266266e+02,
    -8.47708996e+02, 2.98678612e+03,  -9.09210283e+02, -3.46600571e+02, -1.16222662e+04, -5.20351384e+03,
    3.20695153e+04,  -3.84814924e+03, -1.13945860e+04, 1.00000000e+00,  1.00000000e+00,  -5.80933579e+02,
    -6.58619047e+02, 2.10698659e+03,  -2.99608774e+02, -5.66825190e+02, -3.33874282e+03, -2.46326851e+03,
    9.42183488e+03,  -1.64935945e+03, -1.96946409e+03, -4.94493522e+02, -3.77400289e+02, 1.66474610e+03,
    -4.50909788e+02, -3.40942506e+02, -1.02214874e+03, -4.04186409e+02, 2.65520632e+03,  -2.15617134e+02,
    -1.01225404e+03, -1.38578668e+03, -1.20430273e+03, 5.02879092e+03,  -1.38159667e+03, -1.05610484e+03,
    -1.10285014e+03, -1.24230749e+02, 2.58259441e+03,  -2.89052876e+02, -1.06546064e+03, -2.29855672e+02,
    -2.86926088e+02, 1.04169331e+03,  -2.89901479e+02, -2.34010070e+02, -6.35413960e+03, -1.39513122e+03,
    1.49813893e+04,  -2.93805144e+03, -4.29306702e+03, 1.00000000e+00,  1.00000000e+00,  -3.83451064e+02,
    -7.80647125e+02, 2.11993862e+03,  -2.60080932e+02, -6.94759503e+02, -3.23577577e+02, -5.19320285e+02,
    1.69669605e+03,  -2.86939537e+02, -5.65858651e+02, -1.32740169e+02, -1.55585161e+02, 6.57854426e+02,
    -2.43206083e+02, -1.25323014e+02, -7.90604842e+02, -4.37335982e+02, 3.48587312e+03,  -4.76617900e+02,
    -1.78031440e+03, -2.03549076e+02, -6.56923658e+02, 2.05025986e+03,  -6.26845414e+02, -5.61941716e+02,
    -1.41591083e+03, -3.88637468e+02, 3.85761428e+03,  -7.72239044e+02, -1.27982694e+03, -3.97050122e+02,
    -9.78716234e+02, 3.36671040e+03,  -8.72417987e+02, -1.11752606e+03, -7.49998243e+02, -4.32285857e+02,
    1.96443421e+03,  -8.17149821e+01, -6.99435130e+02, 1.00000000e+00,  1.00000000e+00,  1.00000000e+00,
    1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,  1.00000000e+00,
    1.00000000e+00,  1.00000000e+00};
