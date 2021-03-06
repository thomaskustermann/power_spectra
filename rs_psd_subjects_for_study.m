% Execute all code blocks containing strings with subjects and then execute the code for the split desired at bottom

% - Script containts hard-coded divide into different groups of patients (pos vs. neg outcome, split into training and test for prediction and
%   power-spectra statistical computations
% - Training and test sets can be asymmetric as original split was created for connectivity paper and excluded all patients with muscle
%   artifacts (not needed for power-spectra based predictions)
%%
clear
%% patients with comorbidities
subj_com = {
    'p155\d1\process\resting_state\'
    'p156\d1\process\resting_state\'
    'p158\d1\process\resting_state\'
    'p162\d1\process\resting_state\'
    'p210\d1\process\resting_state\'
    'p217\d1\process\resting_state\'
    'p222\d1\process\resting_state\'
    'p225\d1\process\resting_state\'
    'p232\d1\process\resting_state\'
    };
% survivors - overall d1 and training d1 ---> remainder is d1 test survivors
subj_pos_d1 = { % group of patients with positive outcome
    'p151\d1\process\resting_state\'
    'p152\d1\process\resting_state\'
    'p153\d1\process\resting_state\'
    'p157\d1\process\resting_state\'
    'p159\d1\process\resting_state\'
    'p160\d1\process\resting_state\'
    'p163\d1\process\resting_state\'
    'p164\d1\process\resting_state\'
    'p165\d1\process\resting_state\'
    'p168\d1\process\resting_state\'
    'p170\d1\process\resting_state\'
    'p172\d1\process\resting_state\'
    'p173\d1\process\resting_state\'
    'p174\d1\process\resting_state\'
    'p176\d1\process\resting_state\'
    'p177\d1\process\resting_state\'
    'p180\d1\process\resting_state\'
    'p181\d1\process\resting_state\'
    'p183\d1\process\resting_state\'
    'p184\d1\process\resting_state\'
    'p186\d1\process\resting_state\'
    'p187\d1\process\resting_state\'
    'p188\d1\process\resting_state\'
    'p189\d1\process\resting_state\'
    'p190\d1\process\resting_state\'
    'p191\d1\process\resting_state\'
    'p197\d1\process\resting_state\'
    'p198\d1\process\resting_state\'
    'p203\d1\process\resting_state\'
    'p204\d1\process\resting_state\'
    'p207\d1\process\resting_state\'
    'p208\d1\process\resting_state\'
    'p210\d1\process\resting_state\'
    'p215\d1\process\resting_state\'
    'p216\d1\process\resting_state\'
    'p218\d1\process\resting_state\'
    'p219\d1\process\resting_state\'
    'p220\d1\process\resting_state\'
    'p223\d1\process\resting_state\'
    'p224\d1\process\resting_state\'
    'p226\d1\process\resting_state\'
    'p228\d1\process\resting_state\'
    'p229\d1\process\resting_state\'
    'p231\d1\process\resting_state\'
    'p233\d1\process\resting_state\'
    'p236\d1\process\resting_state\'
    'p237\d1\process\resting_state\'
    'p239\d1\process\resting_state\'
    'p241\d1\process\resting_state\'
    'p242\d1\process\resting_state\'
    'p243\d1\process\resting_state\'
    'p244\d1\process\resting_state\'
    'p245\d1\process\resting_state\'
    'p246\d1\process\resting_state\'
    'p247\d1\process\resting_state\'
    'p250\d1\process\resting_state\'
    'p251\d1\process\resting_state\'
    'p254\d1\process\resting_state\'
    'p255\d1\process\resting_state\'
    'p256\d1\process\resting_state\'
    'p260\d1\process\resting_state\'
    'p261\d1\process\resting_state\'
    'p262\d1\process\resting_state\'
    'p263\d1\process\resting_state\'
    'p264\d1\process\resting_state\'
    'p266\d1\process\resting_state\'
    'p268\d1\process\resting_state\'
    'p271\d1\process\resting_state\'
    'p272\d1\process\resting_state\'
    'p275\d1\process\resting_state\'
    'p277\d1\process\resting_state\'
    'p278\d1\process\resting_state\'
    'p279\d1\process\resting_state\'
    'p280\d1\process\resting_state\'
    'p281\d1\process\resting_state\'
    'p282\d1\process\resting_state\'
    'p283\d1\process\resting_state\'
    'p285\d1\process\resting_state\'
    'p289\d1\process\resting_state\'
    'p290\d1\process\resting_state\'
    };
training_pos_d1 = { % training group of patients with positive outcome
    'p159\d1\process\resting_state\'
    'p160\d1\process\resting_state\'
    'p165\d1\process\resting_state\'
    'p172\d1\process\resting_state\'
    'p173\d1\process\resting_state\'
    'p181\d1\process\resting_state\'
    'p184\d1\process\resting_state\'
    'p187\d1\process\resting_state\'
    'p190\d1\process\resting_state\'
    'p191\d1\process\resting_state\'
    'p207\d1\process\resting_state\'
    'p208\d1\process\resting_state\'
    'p216\d1\process\resting_state\'
    'p219\d1\process\resting_state\'
    'p231\d1\process\resting_state\'
    'p239\d1\process\resting_state\'
    'p241\d1\process\resting_state\'
    'p242\d1\process\resting_state\'
    'p243\d1\process\resting_state\'
    'p254\d1\process\resting_state\'
    'p256\d1\process\resting_state\'
    'p261\d1\process\resting_state\'
    'p262\d1\process\resting_state\'
    'p264\d1\process\resting_state\'
    'p277\d1\process\resting_state\'
    'p278\d1\process\resting_state\'
    'p280\d1\process\resting_state\'
    };
% non-survivors - overall d1 and training d1 ---> remainder is d1 test non-survivors
subj_neg_d1 = { % group of patients with negative outcome
    'p154\d1\process\resting_state\'
    'p155\d1\process\resting_state\'
    'p156\d1\process\resting_state\'
    'p158\d1\process\resting_state\'
    'p161\d1\process\resting_state\'
    'p162\d1\process\resting_state\'
    'p166\d1\process\resting_state\'
    'p167\d1\process\resting_state\'
    'p169\d1\process\resting_state\'
    'p171\d1\process\resting_state\'
    'p175\d1\process\resting_state\'
    'p178\d1\process\resting_state\'
    'p182\d1\process\resting_state\'
    'p185\d1\process\resting_state\'
    'p192\d1\process\resting_state\'
    'p193\d1\process\resting_state\'
    'p194\d1\process\resting_state\'
    'p195\d1\process\resting_state\'
    'p196\d1\process\resting_state\'
    'p199\d1\process\resting_state\'
    'p200\d1\process\resting_state\'
    'p201\d1\process\resting_state\'
    'p202\d1\process\resting_state\'
    'p205\d1\process\resting_state\'
    'p206\d1\process\resting_state\'
    'p209\d1\process\resting_state\'
    'p211\d1\process\resting_state\'
    'p212\d1\process\resting_state\'
    'p213\d1\process\resting_state\'
    'p214\d1\process\resting_state\'
    'p217\d1\process\resting_state\'
    'p221\d1\process\resting_state\'
    'p222\d1\process\resting_state\'
    'p225\d1\process\resting_state\'
    'p227\d1\process\resting_state\'
    'p230\d1\process\resting_state\'
    'p232\d1\process\resting_state\'
    'p234\d1\process\resting_state\'
    'p235\d1\process\resting_state\'
    'p238\d1\process\resting_state\'
    'p240\d1\process\resting_state\'
    'p248\d1\process\resting_state\'
    'p249\d1\process\resting_state\'
    'p252\d1\process\resting_state\'
    'p253\d1\process\resting_state\'
    'p257\d1\process\resting_state\'
    'p258\d1\process\resting_state\'
    'p259\d1\process\resting_state\'
    'p265\d1\process\resting_state\'
    'p267\d1\process\resting_state\'
    'p269\d1\process\resting_state\'
    'p270\d1\process\resting_state\'
    'p273\d1\process\resting_state\'
    'p274\d1\process\resting_state\'
    'p276\d1\process\resting_state\'
    'p284\d1\process\resting_state\'
    'p286\d1\process\resting_state\'
    'p287\d1\process\resting_state\'
    };
training_neg_d1 = { % training group of patients with negative outcome
    'p167\d1\process\resting_state\'
    'p169\d1\process\resting_state\'
    'p175\d1\process\resting_state\'
    'p178\d1\process\resting_state\'
    'p193\d1\process\resting_state\'
    'p195\d1\process\resting_state\'
    'p196\d1\process\resting_state\'
    'p201\d1\process\resting_state\'
    'p205\d1\process\resting_state\'
    'p206\d1\process\resting_state\'
    'p209\d1\process\resting_state\'
    'p211\d1\process\resting_state\'
    'p212\d1\process\resting_state\'
    'p235\d1\process\resting_state\'
    'p240\d1\process\resting_state\'
    'p253\d1\process\resting_state\'
    'p265\d1\process\resting_state\'
    'p274\d1\process\resting_state\'
    'p276\d1\process\resting_state\'
    };
% patients in test sample (excluding the ones with muscle or other EEG artifacts)
test_pos_d1_no_muscle = { % patients with positive outcome that show no muscle artifacts
    'p152\d1\process\resting_state\'
    'p153\d1\process\resting_state\'
    'p170\d1\process\resting_state\'
    'p174\d1\process\resting_state\'
    'p177\d1\process\resting_state\'
    'p183\d1\process\resting_state\'
    'p189\d1\process\resting_state\'
    'p197\d1\process\resting_state\'
    'p203\d1\process\resting_state\'
    'p204\d1\process\resting_state\'
    'p215\d1\process\resting_state\'
    'p224\d1\process\resting_state\'
    'p226\d1\process\resting_state\'
    'p228\d1\process\resting_state\'
    'p229\d1\process\resting_state\'
    'p236\d1\process\resting_state\'
    'p245\d1\process\resting_state\'
    'p246\d1\process\resting_state\'
    'p250\d1\process\resting_state\'
    'p251\d1\process\resting_state\'
    'p260\d1\process\resting_state\'
    'p266\d1\process\resting_state\'
    'p268\d1\process\resting_state\'
    'p275\d1\process\resting_state\'
    'p279\d1\process\resting_state\'
    'p281\d1\process\resting_state\'
    'p285\d1\process\resting_state\'
    'p289\d1\process\resting_state\'
    };

test_neg_d1_no_muscle = { % patients with negative outcome that show no muscle artifacts
    'p154\d1\process\resting_state\'
    'p161\d1\process\resting_state\'
    'p171\d1\process\resting_state\'
    'p182\d1\process\resting_state\'
    'p192\d1\process\resting_state\'
    'p194\d1\process\resting_state\'
    'p199\d1\process\resting_state\'
    'p200\d1\process\resting_state\'
    'p213\d1\process\resting_state\'
    'p214\d1\process\resting_state\'
    'p221\d1\process\resting_state\'
    'p238\d1\process\resting_state\'
    'p252\d1\process\resting_state\'
    'p259\d1\process\resting_state\'
    'p267\d1\process\resting_state\'
    'p269\d1\process\resting_state\'
    'p270\d1\process\resting_state\'
    'p284\d1\process\resting_state\'
    };

d1_muscle = { % patients that DO show muscle artifacts
    'p151\d1\process\resting_state\'
    'p155\d1\process\resting_state\'
    'p157\d1\process\resting_state\'
    'p163\d1\process\resting_state\'
    'p164\d1\process\resting_state\'
    'p166\d1\process\resting_state\'
    'p168\d1\process\resting_state\'
    'p176\d1\process\resting_state\'
    'p180\d1\process\resting_state\'
    'p185\d1\process\resting_state\'
    'p186\d1\process\resting_state\'
    'p188\d1\process\resting_state\'
    'p198\d1\process\resting_state\'
    'p202\d1\process\resting_state\'
    'p217\d1\process\resting_state\'
    'p218\d1\process\resting_state\'
    'p220\d1\process\resting_state\'
    'p222\d1\process\resting_state\'
    'p223\d1\process\resting_state\'
    'p227\d1\process\resting_state\'
    'p230\d1\process\resting_state\'
    'p232\d1\process\resting_state\'
    'p233\d1\process\resting_state\'
    'p234\d1\process\resting_state\'
    'p237\d1\process\resting_state\'
    'p244\d1\process\resting_state\'
    'p247\d1\process\resting_state\'
    'p248\d1\process\resting_state\'
    'p249\d1\process\resting_state\'
    'p255\d1\process\resting_state\'
    'p257\d1\process\resting_state\'
    'p258\d1\process\resting_state\'
    'p271\d1\process\resting_state\'
    'p272\d1\process\resting_state\'
    'p282\d1\process\resting_state\'
    'p283\d1\process\resting_state\'
    'p286\d1\process\resting_state\'
    'p287\d1\process\resting_state\'
    'p290\d1\process\resting_state\'
    
    'p273\d1\process\resting_state\' % categorized as "other artifact"
    };


%%%%%%% RUN one of these below (don't execute multiple blocks below) %%%%%%%%
%% split into training and test groups (comorbidities excluded)
test_pos_d1 = subj_pos_d1(~ismember(subj_pos_d1,[subj_com;training_pos_d1]));
test_neg_d1 = subj_neg_d1(~ismember(subj_neg_d1,[subj_com;training_neg_d1]));
%% split into training, test and previously excluded groups (e.g. muscle artifacts; but still excluding comorbidities)
excl_pos_d1 = subj_pos_d1(~ismember(subj_pos_d1,[subj_com;training_pos_d1;test_pos_d1_no_muscle]));
excl_neg_d1 = subj_neg_d1(~ismember(subj_neg_d1,[subj_com;training_neg_d1;test_neg_d1_no_muscle]));
test_pos_d1 = test_pos_d1_no_muscle;
test_neg_d1 = test_neg_d1_no_muscle;
%% split into training and test groups w/o excluding comorbidities or muscle artifacts
test_pos_d1 = subj_pos_d1(~ismember(subj_pos_d1,[training_pos_d1]));
test_neg_d1 = subj_neg_d1(~ismember(subj_neg_d1,[training_neg_d1]));
%% split into training and test groups, including comorbidities
% test_pos_d1 = [test_pos_d1_no_muscle;subj_com(ismember(subj_com,subj_pos_d1))];
% test_neg_d1 = [test_neg_d1_no_muscle;subj_com(ismember(subj_com,subj_neg_d1))];
subj_pos_d1 = subj_pos_d1(~ismember(subj_pos_d1,d1_muscle));
subj_neg_d1 = subj_neg_d1(~ismember(subj_neg_d1,d1_muscle));
test_pos_d1 = subj_pos_d1(~ismember(subj_pos_d1,[d1_muscle;training_pos_d1]));
test_neg_d1 = subj_neg_d1(~ismember(subj_neg_d1,[d1_muscle;training_neg_d1]));
warning('Comorbidity - included; Muscle - excluded')
