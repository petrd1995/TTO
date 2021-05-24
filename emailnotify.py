import smtplib
import os
from dotenv import load_dotenv

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

        subject = 'Vypocet dokoncen!'
        body = 'Vypocet dokoncen!!'

        msg = f'Subject: {subject} \n\n {body}'
        smtp.sendmail(EMAIL_ADDRESS, EMAIL_RECIPIENT, msg)
