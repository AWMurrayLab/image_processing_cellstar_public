import numpy as np
import csv
import pandas as pd
import matplotlib
import os
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import seaborn as sns
import scipy
from scipy import stats
# sns.set(font_scale=2.0)

# this is what we use to name each plot
strain_db={'yFB22':r'pGAL1-WHI5, $\Delta$far1', 'yNTI136':r'pWHI5-WHI5, $\Delta$far1',
           'yFB29':r'$P_{GAL1}-WHI5$', 'yFB30':r'$P_{GAL1}-WHI5$, $\Delta bck2$',
           'yFB41':r'pWHI5-WHI5', 'yFB43':r'$P_{WHI5}-WHI5$', 'yFB25':r'WT', 'yFB86':r'$P_{GAL1}-CLN3$',
          'yFB45':r'pWHI5-WHI5, $\Delta$bck2', 'yFB46':r'$P_{WHI5}-WHI5$, $\Delta bck2$',
           'yFB78':r'pGAL1-WHI5-mVenNB, pACT1-mCherry','yFB79':r'pWHI5-WHI5-mVenNB, pACT1-mCherry',
          'yFB93':r'WT', 'yFB94':r'WT MAT$\alpha$', 'yFB95':r'WT MAT$\alpha$ -leu',
          'yFB96':r'$\Delta$whi5', 'yFB97':r'$\Delta$cln3', 'yFB98':'$\Delta$bck2',
          'yFB99':r'$\Delta$swe1', 'yFB100':r'$\Delta$cln3, $\Delta$whi5',
           'yFB101':r'$\Delta$whi5, $\Delta$bck2',
           'yFB102':r'$\Delta$whi5, $\Delta$cln3, $\Delta$bck2',
           'yFB103':r'$\Delta$whi5, $\Delta$bck2, $\Delta$swe1',
           'yFB104':r'$\Delta$whi5, $\Delta$cln3, $\Delta$bck2, $\Delta$swe1',
          'yFB108':r'$\Delta$whi5, $\Delta$cln3, $\Delta$swe1'}
names = {'mean':'Mean (fL)', 'std':'Standard Deviation (fL)', 'cv':'CV', 'mode':'Mode Volume (fL)',
         'med':'Median Volume (fL)'}

def stats_calc(temp_nums, temp_temp_df, temp_vols, temp_params, temp_ind):
    # adds data to a dataframe based on the experiment in question
    mean_vol = np.sum(temp_nums * temp_vols) / np.sum(temp_nums)
    std_dev = np.sqrt(np.sum(temp_nums * (temp_vols - mean_vol) ** 2) / np.sum(temp_nums))
    cv = std_dev / mean_vol
    #     print np.argmax(temp_nums)
    mode_vol = temp_vols[np.argmax(temp_nums)]

    # calculating the median
    temp_ind1 = np.nonzero(np.cumsum(temp_nums) > np.sum(temp_nums) * 1.0 / 2)[0][0]
    median_vol = temp_vols[temp_ind1]
    temp_vals = [mean_vol, std_dev, cv, mode_vol, median_vol] + temp_params
    # print temp_vals
    # print temp_temp_df.head(10)
    temp_temp_df.loc[temp_ind] = temp_vals
    return temp_temp_df


def stats_calc_v2(temp_nums, temp_temp_df, temp_params, temp_ind):
    # adds data to a dataframe based on the experiment in question
    mean_vol = np.mean(temp_nums)
    std_dev = np.std(temp_nums)
    cv = std_dev / mean_vol
    mode_vol = scipy.stats.mode(temp_nums)
    # calculating the median
    median_vol = np.median(temp_nums)
    temp_vals = [mean_vol, std_dev, cv, mode_vol, median_vol] + temp_params
    temp_temp_df.loc[temp_ind] = temp_vals
    return temp_temp_df


def import_and_plot(temp_df, params, selected_strains, good_copies=False, gal=True):
    # takes as input an empty dataframe with column names assigned and a set of experimental parameters
    # that point to the location of csv files with size distribution data
    # Returns a dataframe filled with statistics about the size distributions and plots of the CDF for the
    # respective datasets
    fig=plt.figure(figsize=[10,6])
    plt.xlabel(r'Volume ($\mu m^3$)')
    plt.ylabel('CDF')

    for i0 in range(len(params)):
        date = params[i0][0]
        strain = params[i0][1]
        csource = params[i0][2]
        gconc = params[i0][3]
        if gal:
            cond = strain+'_'+gconc
        else:
            cond=strain
        sample_num = params[i0][4]
        expt_id = params[i0][5]
        data = []
        with open(bp+date+'_multisizer/'+cond+'/'+expt_id+'.#M3.CSV') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter =',')
            data_count = 0
            line_count = 0
            for row in csv_reader:
                if data_count:
                    line_count+=1
                    data.append(row)
                if len(row)!=0 and row[1]=='um':
                    data_count += 1
        # print data[-2]
        data1 = np.zeros([len(data)-2,4])
        for i1 in range(len(data)-2):
        #     print data[i0]
            data1[i1,:] = np.asarray(data[i1])
        vols = data1[:,1]**3*math.pi/6
        expt_params = [date, cond, csource, gconc, strain]
        cutoff = 2.0  # diameter below which we will not consider the signal as being relevant
        ucutoff = 1000.0  # volume above which we will not consider the signal as being relevant
        ind = np.nonzero(data1[:,1]>cutoff)[0][0]
        ind1 = np.nonzero(vols>ucutoff)[0][0]
#         print ind1, len(vols)
        # plotting
        if i0 in selected_strains:
            if good_copies:
                plt.semilogx(vols[ind:ind1], np.cumsum(data1[ind:ind1,2]/sum(data1[ind:,2])),'-',
                             label = gconc+', '+strain_db[strain])
            else:
                plt.semilogx(vols[ind:ind1], np.cumsum(data1[ind:ind1,2]/sum(data1[ind:,2])),'-', label = cond)
        temp_df = stats_calc(data1[ind:,2], temp_df, vols[ind:], expt_params, i0)
    plt.title('CDF of population size distributions ')
    plt.legend(loc='center left', bbox_to_anchor=(0.7, 0.5))
    return fig, temp_df


def import_and_plot_v4(temp_df, params, selected_strains, temp_path, good_copies=False, gal=True, fs=1.2, label='',
                       lcutoff=10, percent_ucutoff = 0.975):
    # takes as input an empty dataframe with column names assigned and a set of experimental parameters
    # that point to the location of csv files with size distribution data
    # Returns a dataframe filled with statistics about the size distributions
    # cutoff = 2.0  # diameter below which we will not consider the signal as being relevant
    figs=[]
    for i0 in selected_strains:
        date = params[i0][0]
        strain = params[i0][1]
        csource = params[i0][2]
        gconc = params[i0][3]
        if gal:
            cond = strain + '_' + gconc
        else:
            cond = strain
        sample_num = params[i0][4]
        expt_id = params[i0][5]
        data = []
        # print cond
        with open(temp_path + date + '_multisizer/' + cond + '/' + expt_id + '.#M3.CSV') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            data_count = 0
            line_count = 0
            for row in csv_reader:
                if data_count:
                    line_count += 1
                    data.append(row)
                if len(row) != 0 and row[1] == 'um':
                    data_count += 1
        # print data[-2]
        data1 = np.zeros([len(data) - 2, 4])
        for i1 in range(len(data) - 2):
            #     print data[i0]
            data1[i1, :] = np.asarray(data[i1])
        vols = data1[:, 1] ** 3 * math.pi / 6

        temp_pulses = []
        diameters = np.append(data1[:, 1], [data[-2][1]]).astype(float)  # range of diameters measured
        #     print np.diff(diameters)
        for ind2 in range(data1.shape[0]):
            for ind1 in range(int(data[ind2][2])):
                temp_pulses.append((diameters[ind2] + diameters[ind2 + 1]) / 2.0)
        temp_pulses = np.asarray(temp_pulses) ** 3 * math.pi / 6  # pulses in cell volume
        tempval = np.sort(temp_pulses)
        temp_ucutoff = tempval[int(percent_ucutoff*len(tempval))-1]  # determining the upper cutoff relative to the
        # distribution of values for this celltype
        temp_pulses = temp_pulses[np.nonzero(temp_pulses > lcutoff)]
        temp_pulses = temp_pulses[np.nonzero(temp_pulses < temp_ucutoff)]
        expt_params = [date, cond, csource, gconc, strain]
        temp_df = stats_calc_v2(temp_pulses, temp_df, expt_params, i0)
        # plotting
        sns.set(font_scale=2)
        fig=plt.figure(figsize=[12,8])
        sns.distplot(temp_pulses,label=label)
        figs.append(fig)
    return temp_df, figs


def quality_control(temp_df, params, selected_strains, temp_path, good_copies=False, gal=True, fs=1.2, title='', cutoff=2.5):
    # takes as input an empty dataframe with column names assigned and a set of experimental parameters
    # that point to the location of csv files with size distribution data
    # Returns a dataframe filled with statistics about the size distributions
    # cutoff = 2.0  # diameter below which we will not consider the signal as being relevant
    for i0 in selected_strains:
        date = params[i0][0]
        strain = params[i0][1]
        csource = params[i0][2]
        gconc = params[i0][3]
        if gal:
            cond = strain+'_'+gconc
        else:
            cond=strain
        sample_num = params[i0][4]
        expt_id = params[i0][5]
        data = []
        # print cond
        with open(temp_path+date+'_multisizer/'+cond+'/'+expt_id+'.#M3.CSV') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter =',')
            data_count = 0
            line_count = 0
            for row in csv_reader:
                if data_count:
                    line_count+=1
                    data.append(row)
                if len(row)!=0 and row[1]=='um':
                    data_count += 1
        # print data[-2]
        data1 = np.zeros([len(data)-2,4])
        for i1 in range(len(data)-2):
        #     print data[i0]
            data1[i1,:] = np.asarray(data[i1])
        stepsize = data1[1, 1] - data1[0, 1]  # diameter stepsize
        vols = data1[:,1]**3*math.pi/6
        expt_params = [date, cond, csource, gconc, strain]
        # cutoff = 2.0  # diameter below which we will not consider the signal as being relevant
        ucutoff = 1000.0  # volume above which we will not consider the signal as being relevant
        ind = np.nonzero(data1[:,1]>cutoff)[0][0]
        ind1 = np.nonzero(vols>ucutoff)[0][0]
        fig=plt.figure(figsize=[5,5])
        plt.semilogx(vols[ind:ind1], np.cumsum(data1[ind:ind1,2]/(sum(data1[ind:ind1,2]))),
                     '-',label = gconc+', '+strain_db[strain])
        plt.title(date+' '+strain+' '+csource+' '+gconc)
        fig.savefig(temp_path+date+'_multisizer/'+cond+'/'+cond+'{0}_CDF.png'.format(params[i0][4]))
        fig=plt.figure(figsize=[5,5])
        plt.semilogx(vols[ind:ind1], data1[ind:ind1,2]/(sum(data1[ind:ind1,2])*stepsize*math.pi*data1[ind:ind1,1]**2/2),
                     '-',label = gconc+', '+strain_db[strain])
        plt.title(date+' '+strain+' '+csource+' '+gconc)
        fig.savefig(temp_path+date+'_multisizer/'+cond+'/'+cond+'{0}_PDF.png'.format(params[i0][4]))
        # plotting
        temp_df = stats_calc(data1[ind:,2], temp_df, vols[ind:], expt_params, i0)
    return temp_df


def import_boxcox(temp_df, params, selected_strains, temp_path, good_copies=False, gal=True, fs=1.2, title=''):
    sns.set(font_scale=fs)
    # takes as input an empty dataframe with column names assigned and a set of experimental parameters
    # that point to the location of csv files with size distribution data
    # Returns a dataframe filled with statistics about the size distributions and plots of the CDF for the
    # respective datasets
    fig=plt.figure(figsize=[10,6])
    plt.xlabel(r'Box Cox Volume ($\mu m^3$)')
    plt.ylabel('PDF')
    temp_vals=[]
    for i0 in selected_strains:
        date = params[i0][0]
        strain = params[i0][1]
        csource = params[i0][2]
        gconc = params[i0][3]
        if gal:
            cond = strain+'_'+gconc
        else:
            cond=strain
        sample_num = params[i0][4]
        expt_id = params[i0][5]
        data = []
        with open(temp_path+date+'_multisizer/'+cond+'/'+expt_id+'.#M3.CSV') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter =',')
            data_count = 0
            line_count = 0
            for row in csv_reader:
                if data_count:
                    line_count+=1
                    data.append(row)
                if len(row)!=0 and row[1]=='um':
                    data_count += 1
        # print data[-2]
        data1 = np.zeros([len(data)-2,4])
        for i1 in range(len(data)-2):
        #     print data[i0]
            data1[i1,:] = np.asarray(data[i1])
        vols = data1[:,1]**3*math.pi/6
        expt_params = [date, cond, csource, gconc, strain]
        cutoff = 2.0  # diameter below which we will not consider the signal as being relevant
        ucutoff = 1000.0  # volume above which we will not consider the signal as being relevant
        ind = np.nonzero(data1[:,1]>cutoff)[0][0]
        ind1 = np.nonzero(vols>ucutoff)[0][0]
        temp1 = 0
        vals =np.array([])
        for temp1 in range(ind, ind1):
            vals = np.append(vals, vols[temp1]*np.ones([int(data1[temp1,2])]))
#         print ind1, len(vols)
        # plotting
        vals1=scipy.stats.boxcox(vals)
        # print vals1.shape
        if i0 in selected_strains:
            if good_copies:
                sns.distplot(vals1[0],label = gconc+', '+strain_db[strain])
            else:
                sns.distplot(vals1[0], label=gconc + ', ' + strain_db[strain])
        temp_df = stats_calc(data1[ind:,2], temp_df, vols[ind:], expt_params, i0)
        temp_vals.append(vals1)
    plt.title('PDF of BOXCOX size distributions '+title)
    plt.legend(loc='center left', bbox_to_anchor=(0.6, 0.5))
    return fig, temp_df, temp_vals


def import_beads(temp_df, params, selected_strains, temp_path, good_copies=False, gal=True, fs=1.2, title='', cutoff=1.0):
    # takes as input an empty dataframe with column names assigned and a set of experimental parameters
    # that point to the location of csv files with size distribution data
    # Returns a dataframe filled with statistics about the size distributions
    # cutoff = 2.0  # diameter below which we will not consider the signal as being relevant
    print selected_strains
    fig = plt.figure(figsize=[10,6])
    for i0 in selected_strains:
        date = params[i0][0]
        strain = params[i0][1]
        csource = params[i0][2]
        gconc = params[i0][3]
        cond=strain
        sample_num = params[i0][4]
        expt_id = params[i0][5]
        data = []
        # print cond
        with open(temp_path+date+'_multisizer/'+strain+'/'+expt_id+'.#M3.CSV') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter =',')
            data_count = 0
            line_count = 0
            for row in csv_reader:
                if data_count:
                    line_count+=1
                    data.append(row)
                if len(row)!=0 and row[1]=='um':
                    data_count += 1
        data1 = np.zeros([len(data)-2,4])
        for i1 in range(len(data)-2):
            data1[i1,:] = np.asarray(data[i1])
        vols = data1[:,1]**3*math.pi/6
        ind = np.nonzero(data1[:, 1] > cutoff)[0][0]
        pdf = data1[ind:,2]/(np.sum(data1[ind:,2])*math.pi*data1[ind:,1]**2/2)
        plt.semilogx(vols[ind:], pdf, '-', label=date+', '+ strain+ ', ' + str(sample_num))
        plt.title('PDF of population size distributions ' + title)
        plt.xlim(xmin=10.0, xmax=30.0)
        plt.legend()
        expt_params = [date, cond, csource, gconc, strain]
        temp_df = stats_calc(data1[ind:, 2], temp_df, vols[ind:], expt_params, i0)
    return fig, temp_df


def import_data(params, selected_strains, temp_path, gal=True, lcutoff=10, percent_ucutoff = 0.975):
    # takes as input an empty dataframe with column names assigned and a set of experimental parameters
    # that point to the location of csv files with size distribution data
    # Returns a dataframe filled with the various data points. This dataframe will be v large.
    temp_cols = ['Strain', 'Genotype', 'Galactose Concentration', 'Volume', 'Sample Number']
    temp_df=pd.DataFrame(columns=temp_cols)
    for i0 in selected_strains:
        date = params[i0][0]
        strain = params[i0][1]
        csource = params[i0][2]
        gconc = params[i0][3]
        if gal:
            cond = strain + '_' + gconc
        else:
            cond = strain
        sample_num = params[i0][4]
        expt_id = params[i0][5]
        data = []
        # print cond
        with open(temp_path + date + '_multisizer/' + cond + '/' + expt_id + '.#M3.CSV') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            data_count = 0
            line_count = 0
            for row in csv_reader:
                if data_count:
                    line_count += 1
                    data.append(row)
                if len(row) != 0 and row[1] == 'um':
                    data_count += 1
        # print data[-2]
        data1 = np.zeros([len(data) - 2, 4])
        for i1 in range(len(data) - 2):
            #     print data[i0]
            data1[i1, :] = np.asarray(data[i1])
        temp_pulses = []
        diameters = np.append(data1[:, 1], [data[-2][1]]).astype(float)  # range of diameters measured
        for ind2 in range(data1.shape[0]):
            for ind1 in range(int(data[ind2][2])):
                temp_pulses.append((diameters[ind2] + diameters[ind2 + 1]) / 2.0)
        temp_pulses = np.asarray(temp_pulses) ** 3 * math.pi / 6  # pulses in cell volume
        tempval = np.sort(temp_pulses)
        # temp_ucutoff = tempval[int(percent_ucutoff*len(tempval))-1]  # determining the upper cutoff relative to the
        # distribution of values for this celltype
        np.save(temp_path + date + '_multisizer/' + cond + '/' + expt_id,temp_pulses)
        if not os.path.exists(temp_path + date + '_multisizer/csv_pulses/'):
            os.makedirs(temp_path + date + '_multisizer/csv_pulses/')
        out_path = temp_path + date + '_multisizer/csv_pulses/' + expt_id+'.csv'
        np.savetxt(out_path,temp_pulses)
        temp_pulses = temp_pulses[np.nonzero(temp_pulses > lcutoff)]
        # temp_pulses = temp_pulses[np.nonzero(temp_pulses < temp_ucutoff)]

        gconc = [params[i0][3] for ind in range(len(temp_pulses))]
        strain = [params[i0][1] for ind in range(len(temp_pulses))]
        genotype = [strain_db[params[i0][1]] for ind in range(len(temp_pulses))]
        run_num = [i0 for ind in range(len(temp_pulses))]
        temp_data = zip(*[strain,genotype,gconc,temp_pulses, run_num])
        temp_df1 = pd.DataFrame(temp_data, columns=temp_cols)
        temp_df=temp_df.append(temp_df1)
    return temp_df


def import_and_plot_v5(temp_df, params, selected_strains, temp_path, strain_labels, gal=True,
                       lcutoff=8, percent_ucutoff = 0.99):
    # takes as input an empty dataframe with column names assigned and a set of experimental parameters
    # that point to the location of csv files with size distribution data
    # Returns a dataframe filled with statistics about the size distributions
    # cutoff = 2.0  # diameter below which we will not consider the signal as being relevant
    plot_vals = []
    labels=[]
    for i0 in selected_strains:
        date = params[i0][0]
        strain = params[i0][1]
        csource = params[i0][2]
        gconc = params[i0][3]
        if gal:
            cond = strain + '_' + gconc
        else:
            cond = strain
        sample_num = params[i0][4]
        expt_id = params[i0][5]
        data = []
        # print cond
        with open(temp_path + date + '_multisizer/' + cond + '/' + expt_id + '.#M3.CSV') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            data_count = 0
            line_count = 0
            for row in csv_reader:
                if data_count:
                    line_count += 1
                    data.append(row)
                if len(row) != 0 and row[1] == 'um':
                    data_count += 1
        # print data[-2]
        data1 = np.zeros([len(data) - 2, 4])
        for i1 in range(len(data) - 2):
            #     print data[i0]
            data1[i1, :] = np.asarray(data[i1])
        vols = data1[:, 1] ** 3 * math.pi / 6

        temp_pulses = []
        diameters = np.append(data1[:, 1], [data[-2][1]]).astype(float)  # range of diameters measured
        #     print np.diff(diameters)
        for ind2 in range(data1.shape[0]):
            for ind1 in range(int(data[ind2][2])):
                temp_pulses.append((diameters[ind2] + diameters[ind2 + 1]) / 2.0)
        temp_pulses = np.asarray(temp_pulses) ** 3 * math.pi / 6  # pulses in cell volume
        tempval = np.sort(temp_pulses)
        # temp_ucutoff = tempval[int(percent_ucutoff*len(tempval))-1]  # determining the upper cutoff relative to the
        # distribution of values for this celltype
        temp_pulses = temp_pulses[np.nonzero(temp_pulses > lcutoff)]
        # temp_pulses = temp_pulses[np.nonzero(temp_pulses < temp_ucutoff)]
        temp1,temp2=np.unique(temp_pulses,return_counts=True)
        temp2=1.0*np.cumsum(temp2)/np.sum(temp2)
        plot_vals.append([temp1,temp2])
        print len(plot_vals[-1][0]), data1.shape[0]
        labels.append(strain_labels[strain])
        # print np.amax(plot_vals[-1][1])
    # plotting
    sns.set(font_scale=2)
    fig=plt.figure(figsize=[8,5])
    for i0 in range(len(selected_strains)):
        plt.semilogx(plot_vals[i0][0], plot_vals[i0][1],label=labels[i0])
    plt.legend(loc=[0.45,0.1])
    plt.xlabel('Volume (fL)')
    plt.xlim(xmax=200)
    plt.ylabel('CDF')
    return temp_df, fig
