import pandas as pd
import glob


img = ['01', '02', '03', '04', '05']
ker = ['3', '5', '25', '49', '99']
proc = ['4', '8', '16', '32']


#print(kernel_time, file_time, conv_time, save_time)


for i in img:
    for k in ker:
        out_put_file = open("./clean_data/MPI_{}_{}".format(i,k),"w")
        out_put_file.write("TYPE,T.KERNEL,T.FILE,T.CONV,T.SAVE,THREAD\n")
        #SERIAL
        filename = "./SERIAL/results/SERIAL_{}_{}.*".format(i,k)
        filename = glob.glob(filename)[0]
        raw_data = pd.read_csv(filename, names = ["PRO","ACTION","TIME"])
        kernel = raw_data[raw_data['ACTION']==' READKERNEL']
        file = raw_data[raw_data['ACTION']==' READFILE']
        conv = raw_data[raw_data['ACTION']==' CONV']
        save = raw_data[raw_data['ACTION']==' SAVE']
        kernel_time = kernel.TIME.sum()
        file_time = file.TIME.sum()
        conv_time = conv.TIME.sum()
        save_time = save.TIME.sum()
        out_put_file.write("SERIAL,{},{},{},{},{}\n".format(kernel_time,file_time,conv_time,save_time,1))
        #MPI
        for p in proc:
            filename = "./MPI/results/MPI_{}_{}_{}.*".format(i,k,p)
            filename = glob.glob(filename)[0]
            raw_data = pd.read_csv(filename, names = ["PRO","ACTION","TIME"])
            kernel = raw_data[raw_data['ACTION']==' READKERNEL']
            file = raw_data[raw_data['ACTION']==' READFILE']
            conv = raw_data[raw_data['ACTION']==' CONV']
            save = raw_data[raw_data['ACTION']==' SAVE']
            kernel_time = kernel.TIME.sum()
            file_time = file.TIME.sum()
            conv_time = conv.TIME.sum()
            save_time = save.TIME.sum()
            conv_time = conv_time - file_time - save_time
            out_put_file.write("MPI,{},{},{},{},{}\n".format(kernel_time,file_time,conv_time,save_time,p))