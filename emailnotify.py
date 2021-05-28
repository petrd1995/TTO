import smtplib
import os
from dotenv import load_dotenv
'''
    Function used to notify one, that the optimization has ended. 
    Has to be set up by creating .env file in the same folder as the rest of the files and there one has to define EMAIL_NAME = your email, EMAIL_PASSWORD = your email's password and EMAIL_REC = recipient's email (can be just EMAIL_REC = EMAIL_NAME).
    It is set up only for gmail, and in order to access it from the program, you have to allow third party access to your email in your settings.
'''
load_dotenv()

EMAIL_ADDRESS = os.getenv('EMAIL_NAME')
EMAIL_PASS = os.getenv('EMAIL_PASSWORD')
EMAIL_RECIPIENT = os.getenv('EMAIL_REC')

def notify():
    with smtplib.SMTP('smtp.gmail.com', 587) as smtp:
        smtp.ehlo()
        smtp.starttls()
        smtp.ehlo()

        smtp.login(EMAIL_ADDRESS, EMAIL_PASS)

        subject = 'Computation finished!'
        body = 'Computation finished!!'

        msg = f'Subject: {subject} \n\n {body}'
        smtp.sendmail(EMAIL_ADDRESS, EMAIL_RECIPIENT, msg)
