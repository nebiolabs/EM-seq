<<<<<<< HEAD
# Capfile
require "capistrano/setup"
require "capistrano/deploy"
require "capistrano/scm/git"
require 'capistrano/bundler'
require 'capistrano/rbenv'
require 'dotenv'
Dotenv.load
install_plugin Capistrano::SCM::Git

set :rbenv_type, :user
set :rbenv_ruby, '3.4.5'



=======
# Load DSL and set up stages
require "capistrano/setup"

# Include default deployment tasks
require "capistrano/deploy"

require "capistrano/scm/git"
install_plugin Capistrano::SCM::Git


# Load custom tasks from `lib/capistrano/tasks` if you have any defined
Dir.glob("lib/capistrano/tasks/*.rake").each { |r| import r }
>>>>>>> master
