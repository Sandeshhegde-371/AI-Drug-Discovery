import { Link } from 'react-router-dom';
import { BookOpen, FlaskConical, FileSpreadsheet, ChevronRight, Sparkles, Brain, Zap } from 'lucide-react';
import backgroundImage from 'figma:asset/001ea9ce6905b4f7224a4b7530fc75c8991444c0.png';

export function HomePage() {
  return (
    <div className="min-h-screen bg-[#0a0a0f] text-white relative overflow-hidden">
      {/* Animated Background */}
      <div className="absolute inset-0">
        <div 
          className="absolute inset-0 opacity-15"
          style={{
            backgroundImage: `url(${backgroundImage})`,
            backgroundSize: 'cover',
            backgroundPosition: 'center',
            backgroundRepeat: 'no-repeat'
          }}
        />
        {/* Gradient Overlays */}
        <div className="absolute top-0 right-0 w-[600px] h-[600px] bg-purple-600/20 rounded-full blur-[120px]" />
        <div className="absolute bottom-0 left-0 w-[500px] h-[500px] bg-blue-600/15 rounded-full blur-[100px]" />
        <div className="absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 w-[400px] h-[400px] bg-purple-500/10 rounded-full blur-[80px]" />
      </div>
      
      {/* Content */}
      <div className="relative z-10">
        {/* Header */}
        <header className="px-8 py-6 flex items-center justify-between backdrop-blur-sm bg-white/5 border-b border-white/10">
          <div className="flex items-center gap-3">
            <div className="p-2 bg-gradient-to-br from-purple-500 to-purple-700 rounded-lg">
              <Brain className="w-6 h-6 text-white" />
            </div>
            <div className="flex flex-col">
              <span className="text-xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent">
                AI-DRIVEN DDD
              </span>
              <span className="text-xs text-gray-500">Drug Discovery Platform</span>
            </div>
          </div>
          
          <button
  onClick={() => {
    window.location.href = "http://localhost:3001"; // change this to your AI app's local or deployed URL
  }}
  className="px-6 py-2.5 rounded-full bg-gradient-to-r from-purple-500 to-purple-700 text-white hover:shadow-lg hover:shadow-purple-500/50 transition-all hover:scale-105"
>
  Get Started
</button>

        </header>

        {/* Main Content */}
        <main className="container mx-auto px-8 py-16">
          <div className="max-w-7xl mx-auto">
            {/* Hero Section */}
            <div className="text-center space-y-8 mb-20">
              <div className="inline-flex items-center gap-2 px-4 py-2 rounded-full bg-purple-500/10 border border-purple-500/20 mb-4">
                <Sparkles className="w-4 h-4 text-purple-400" />
                <span className="text-sm text-purple-300">Powered by Advanced AI & Machine Learning</span>
              </div>
              
              <h1 className="text-7xl leading-tight">
                <span className="block text-white mb-2">Transform</span>
                <span className="bg-gradient-to-r from-purple-400 via-purple-500 to-blue-500 bg-clip-text text-transparent">
                  Drug Discovery
                </span>
              </h1>
              
              <p className="text-xl text-gray-400 max-w-3xl mx-auto leading-relaxed">
                Leverage cutting-edge AI pipelines to generate new molecules, predict properties, 
                and accelerate your research from concept to candidate
              </p>

              {/* Quick Stats */}
              <div className="flex items-center justify-center gap-12 pt-8">
                <div className="text-center">
                  <div className="text-4xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent mb-1">4</div>
                  <div className="text-sm text-gray-500">AI Pipelines</div>
                </div>
                <div className="w-px h-12 bg-white/10" />
                <div className="text-center">
                  <div className="text-4xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent mb-1">∞</div>
                  <div className="text-sm text-gray-500">Molecules</div>
                </div>
                <div className="w-px h-12 bg-white/10" />
                <div className="text-center">
                  <div className="text-4xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent mb-1">10x</div>
                  <div className="text-sm text-gray-500">Faster</div>
                </div>
              </div>
            </div>

            {/* Navigation Cards */}
            <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
              <Link 
                to="/overview"
                className="group relative bg-gradient-to-br from-purple-950/50 to-purple-900/30 border border-purple-500/20 rounded-3xl p-8 hover:border-purple-400/50 transition-all duration-500 hover:shadow-2xl hover:shadow-purple-500/20 overflow-hidden"
              >
                {/* Hover Glow Effect */}
                <div className="absolute inset-0 bg-gradient-to-br from-purple-500/0 to-purple-500/0 group-hover:from-purple-500/10 group-hover:to-transparent transition-all duration-500" />
                
                <div className="relative z-10">
                  <div className="mb-6 p-4 bg-gradient-to-br from-purple-500 to-purple-700 rounded-2xl w-fit group-hover:scale-110 group-hover:rotate-3 transition-all duration-500">
                    <BookOpen className="w-8 h-8 text-white" />
                  </div>
                  
                  <h2 className="text-3xl text-white mb-3 group-hover:text-purple-300 transition-colors">Overview</h2>
                  <p className="text-gray-400 leading-relaxed mb-6">
                    Discover how our AI-driven platform revolutionizes molecular design and drug discovery workflows
                  </p>
                  
                  <div className="flex items-center gap-2 text-purple-400 group-hover:gap-4 transition-all">
                    <span>Learn More</span>
                    <ChevronRight className="w-4 h-4" />
                  </div>
                </div>
              </Link>

              <Link 
                to="/uses"
                className="group relative bg-gradient-to-br from-blue-950/50 to-blue-900/30 border border-blue-500/20 rounded-3xl p-8 hover:border-blue-400/50 transition-all duration-500 hover:shadow-2xl hover:shadow-blue-500/20 overflow-hidden"
              >
                <div className="absolute inset-0 bg-gradient-to-br from-blue-500/0 to-blue-500/0 group-hover:from-blue-500/10 group-hover:to-transparent transition-all duration-500" />
                
                <div className="relative z-10">
                  <div className="mb-6 p-4 bg-gradient-to-br from-blue-500 to-blue-700 rounded-2xl w-fit group-hover:scale-110 group-hover:rotate-3 transition-all duration-500">
                    <FlaskConical className="w-8 h-8 text-white" />
                  </div>
                  
                  <h2 className="text-3xl text-white mb-3 group-hover:text-blue-300 transition-colors">Uses</h2>
                  <p className="text-gray-400 leading-relaxed mb-6">
                    Explore four powerful AI pipelines for generation, screening, prediction, and prioritization
                  </p>
                  
                  <div className="flex items-center gap-2 text-blue-400 group-hover:gap-4 transition-all">
                    <span>Explore Pipelines</span>
                    <ChevronRight className="w-4 h-4" />
                  </div>
                </div>
              </Link>

              <Link 
                to="/sample"
                className="group relative bg-gradient-to-br from-indigo-950/50 to-indigo-900/30 border border-indigo-500/20 rounded-3xl p-8 hover:border-indigo-400/50 transition-all duration-500 hover:shadow-2xl hover:shadow-indigo-500/20 overflow-hidden"
              >
                <div className="absolute inset-0 bg-gradient-to-br from-indigo-500/0 to-indigo-500/0 group-hover:from-indigo-500/10 group-hover:to-transparent transition-all duration-500" />
                
                <div className="relative z-10">
                  <div className="mb-6 p-4 bg-gradient-to-br from-indigo-500 to-indigo-700 rounded-2xl w-fit group-hover:scale-110 group-hover:rotate-3 transition-all duration-500">
                    <FileSpreadsheet className="w-8 h-8 text-white" />
                  </div>
                  
                  <h2 className="text-3xl text-white mb-3 group-hover:text-indigo-300 transition-colors">Sample</h2>
                  <p className="text-gray-400 leading-relaxed mb-6">
                    View comprehensive sample inputs and outputs for each AI pipeline with real examples
                  </p>
                  
                  <div className="flex items-center gap-2 text-indigo-400 group-hover:gap-4 transition-all">
                    <span>View Samples</span>
                    <ChevronRight className="w-4 h-4" />
                  </div>
                </div>
              </Link>
            </div>

            {/* Features Section */}
            <div className="mt-24 grid grid-cols-1 md:grid-cols-3 gap-8">
              <div className="text-center space-y-4">
                <div className="p-4 bg-purple-500/10 rounded-2xl w-fit mx-auto border border-purple-500/20">
                  <Sparkles className="w-8 h-8 text-purple-400" />
                </div>
                <h3 className="text-xl text-white">Generate Molecules</h3>
                <p className="text-gray-500">Create novel molecular structures using deep generative AI models</p>
              </div>
              
              <div className="text-center space-y-4">
                <div className="p-4 bg-blue-500/10 rounded-2xl w-fit mx-auto border border-blue-500/20">
                  <Brain className="w-8 h-8 text-blue-400" />
                </div>
                <h3 className="text-xl text-white">Predict Properties</h3>
                <p className="text-gray-500">Estimate activity, toxicity, and drug-likeness with ML models</p>
              </div>
              
              <div className="text-center space-y-4">
                <div className="p-4 bg-indigo-500/10 rounded-2xl w-fit mx-auto border border-indigo-500/20">
                  <Zap className="w-8 h-8 text-indigo-400" />
                </div>
                <h3 className="text-xl text-white">Accelerate Research</h3>
                <p className="text-gray-500">Reduce costs while dramatically speeding up discovery</p>
              </div>
            </div>
          </div>
        </main>

        {/* Footer */}
        <footer className="border-t border-white/10 mt-20 py-8 px-8 text-center text-gray-500 text-sm">
          <p>© 2025 AI-Driven Drug Discovery Platform. Accelerating the future of medicine.</p>
        </footer>
      </div>
    </div>
  );
}
