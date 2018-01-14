import smtplib

server = 'pikachu.udel.edu'
db = 'MG_priv_BSseq'
 
def sendmail():
        print('Sending job complete mail')
        to = 'kakrana@gmail.com'
        gmail_user = 'daemon.blake.lab@gmail.com'
        gmail_pwd = 'WeLComE@MeYersLab!'
        smtpserver = smtplib.SMTP("smtp.gmail.com",587)
        smtpserver.ehlo()
        smtpserver.starttls()
        smtpserver.ehlo
        smtpserver.login(gmail_user, gmail_pwd)
        header = 'To:' + to + '\n' + 'From: ' + gmail_user + '\n' + 'Subject:Script run finished \n'
        #print (header)
        msg = (header + 'Master, \nYour BS seq clustering and methylation analysis script has just finished run at:%s on :%s\n\n' % (server,db))
        #msg = (header + '\n This is test msg from your server \n\n')
        smtpserver.sendmail(gmail_user, to, msg)
        print ('Mail Sent!')
        smtpserver.close()      
#####
sendmail()