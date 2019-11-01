pipeline {
    agent { docker { image 'pmantini/assignment-cosc6380:latest' } }

    environment {
        PATH = "env/bin/:$PATH"
    }
    stages {
        stage('build') {
            steps {
                sh 'python dip_hw4_filter.py -i Lenna.png -n gaussian -v 100'
                sh 'python dip_hw4_filter.py -i Lenna.png -n gaussian -s 5 -v 100'
                sh 'python dip_hw4_filter.py -i Lenna.png -n bipolar -npa 0.1 -npb 0 -m max'
                sh 'python dip_hw4_filter.py -i Lenna.png -n bipolar -npa 0 -npb 0.1 -m min'
                sh 'python dip_hw4_filter.py -i Lenna.png -n bipolar -npa 0.1 -npb 0.1'
                sh 'python dip_hw4_filter.py -i Lenna.png -n bipolar -npa 0.1 -npb 0.1 -m alpha_trimmed -s 5 -p 5'
                sh 'python dip_hw4_filter.py -i Lenna.png -n gaussian -v 100 -m arithmetic_mean'
                sh 'python dip_hw4_filter.py -i Lenna.png -n gaussian -v 100 -m geometric_mean'
                sh 'python dip_hw4_filter.py -i Lenna.png -n bipolar -npa 0 -npb 0.1 -m contra_harmonic -o -1.5'
                sh 'python dip_hw4_filter.py -i Lenna.png -n bipolar -npa 0.1 -npb 0 -m contra_harmonic'
            }
        }
    }
    post {
        always {
            archiveArtifacts artifacts: 'output/**/*.* ', onlyIfSuccessful: true
        }
    }
}